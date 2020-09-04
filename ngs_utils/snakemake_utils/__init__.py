import math
from random import random
import string
import os
import subprocess
import sys
from datetime import datetime
from os.path import dirname, abspath, join, isdir
from ngs_utils import logger
from ngs_utils.file_utils import safe_mkdir
from ngs_utils.utils import flatten
import tempfile
import yaml
from ngs_utils.call_process import run_simple


def package_path():
    return dirname(abspath(__file__))


def get_submit_script():
    return 'python ' + join(package_path(), 'submit')


def make_cluster_cmdl(log_dir, refdata, app_name=''):
    """ Generates cluster command line parameters for snakemake
    """
    if not refdata.cluster_cmd:
        logger.critical(f'Automatic cluster submission is not supported for the machine "{refdata.name}"')

    cluster_submitter = get_submit_script()
    timestamp = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    cluster_cmdl = f' --cluster "{cluster_submitter} {timestamp} {log_dir} {app_name}"'

    # Also overriding jobscript?
    jobscript = refdata.cluster_jobscript
    if jobscript:
        safe_mkdir(log_dir)
        jobscript_file = join(log_dir, 'jobscript.sh')
        with open(jobscript_file, 'w') as f_out:
            f_out.write(jobscript.replace('{path}', os.environ["PATH"]))
        cluster_cmdl += f' --jobscript "{jobscript_file}"'

    return cluster_cmdl


DEFAULT_RESTART_TIMES = 2

def run_snakemake(snakefile, conf, cores=None, output_dir=None, forcerun=None,
                  unlock=False, dryrun=False, target_rules=None, debug=False,
                  log_dir=None, dag=None, report=None, restart_times=1,
                  tibanna_cfg=None,
                  resources=None, cluster_param=None, cluster_log_dir=None,
                  local_cores=None, ncpus_per_batch=None, ncpus_per_sample=None,
                  tmp_dirs:list = None):

    ##########################
    #### Preparing config ####
    ##########################

    if unlock: conf['unlock'] = 'yes'

    if debug:
        conf['debug'] = 'yes'
        if restart_times is None:
            restart_times = 0
    if restart_times is None:
        restart_times = DEFAULT_RESTART_TIMES
    restart_times = int(restart_times)

    if ncpus_per_batch:
        conf['threads_per_batch'] = ncpus_per_batch
    if ncpus_per_sample:
        conf['threads_per_sample'] = ncpus_per_sample

    if log_dir:
        safe_mkdir(log_dir)
        conf_f = open(join(log_dir, '.conf.yaml'), 'w')
    else:
        conf_f = tempfile.NamedTemporaryFile(mode='wt', delete=False)
    yaml.dump(conf, conf_f)
    conf_f.close()

    ###############################
    #### Building command line ####
    ###############################

    if forcerun:
        forcerun = " ".join(forcerun.split(','))

    tibanna_opts = ''
    if tibanna_cfg:
        output_s3 = tibanna_cfg['output_s3']
        output_bucket_name = output_s3.split('/')[0]
        if ':' in output_bucket_name:
            output_bucket_name = output_bucket_name.split(':')[1]
        step_func_name = setup_tibanna(tibanna_cfg['id'], [output_bucket_name])
        tibanna_opts = f'--tibanna --default-remote-prefix {output_s3} --tibanna-sfn {step_func_name}'

    cmd = (
        f'snakemake '
        f'{" ".join(flatten([target_rules])) if target_rules else ""} ' +
        f'--snakefile {snakefile} '
        f'--printshellcmds '
        f'{"--dryrun " if dryrun else ""}'
        f'--rerun-incomplete '
        f'{"--dag " if dag else ""}'
        f'{f"--report {report} " if report else ""}'
        f'{f"--directory {output_dir} " if output_dir else ""}'
        f'--cores {cores} '
        f'{f"--local-cores {local_cores} " if local_cores else ""}'
        f'{f"--restart-times {restart_times - 1} " if restart_times > 1 else ""}'
        f'{cluster_param if cluster_param else ""} '
        f'--configfile {conf_f.name} ' +
        f'{f"--forcerun {forcerun} " if forcerun else ""}' +
        f'{f"--resources {resources} " if resources else ""} '
        f'{tibanna_opts}'
    )

    #################
    #### Running ####
    #################

    if unlock:
        print('* Unlocking previous run... *')
        run_simple(cmd + ' --unlock')
        print('* Now rerunning *')

    try:
        run_simple(cmd)
    except subprocess.CalledProcessError:
        logger.error('--------')
        logger.error(f'Error: snakemake returned a non-zero status. Working directory: {output_dir}')
        if cluster_log_dir and isdir(cluster_log_dir):
            run_simple(f'chmod -R a+r {cluster_log_dir}', silent=True)
            logger.error(f'Review cluster job logs in {cluster_log_dir}')
        for tmp_dir in tmp_dirs or []: tmp_dir.cleanup()
        sys.exit(1)
    except KeyboardInterrupt:
        logger.error('--------')
        logger.error(f'Interrupted. Fixing logs permissions. Working directory: {output_dir}')
        if cluster_log_dir and isdir(cluster_log_dir):
            run_simple(f'chmod -R a+r {cluster_log_dir}', silent=True)
            logger.error(f'Review cluster job logs in {cluster_log_dir}')
        for tmp_dir in tmp_dirs or []: tmp_dir.cleanup()
        sys.exit(1)
    else:
        logger.info('--------')
        if cluster_log_dir and isdir(cluster_log_dir):
            run_simple(f'chmod -R a+r {cluster_log_dir}', silent=True)
        logger.info(f'Finished. Output directory: {output_dir}')
        for tmp_dir in tmp_dirs or []: tmp_dir.cleanup()


def prep_resources(num_batches=None, num_samples=None,
                   ncpus_requested=None, is_cluster=False, is_silent=False,
                   ncpus_per_node=None
                   ):
    """ Determines the number of cpus used by a job and the total number of cpus
        available to snakemake scheduler.
        :returns ncpus_per_batch, ncpus_per_sample, ncpus_available, ncpus_per_node=None
    """
    # Checking presets for known HPC clusters, otherwise assuming a for single-machine AWS or local run
    # and just taking the number of available CPUs:

    class Resources:
        def __init__(self):
            self.ncpus_per_batch = None
            self.ncpus_per_sample = None
            self.ncpus_available = None
            self.ncpus_per_node = None
            self.ncpus_local = None

    resources = Resources()

    if is_cluster:
        # we are not resticted to one machine, so can submit many jobs and let the scheduler figure out the queue
        resources.ncpus_available = ncpus_requested or ncpus_per_node
        resources.ncpus_per_node = ncpus_per_node
        if ncpus_per_node:
            logger.info(f'Number of CPUs on a cluster node: {ncpus_per_node}')
        resources.ncpus_local = 1
    else:
        try:
            ncpus_on_this_machine = len(os.sched_getaffinity(0))
        except:
            ncpus_on_this_machine = os.cpu_count()
        if ncpus_on_this_machine:
            logger.info(f'Number of CPUs on this machine : {ncpus_on_this_machine}')
        # scheduling is on Snakemake, so restricting to the number of available cpus on this machine
        resources.ncpus_available = min(ncpus_on_this_machine, ncpus_requested or math.inf)
        resources.ncpus_per_node = None
        resources.ncpus_local = ncpus_on_this_machine

    def adjust_ncpus_per_job(ncpus, max_ncpus_per_job=10, msg=''):
        """ Adjusting the number of cpus to a number below <max_ncpus_per_job>.
            Say, if we have more than 20 cpus on a node and only 1 batch, we should adjust
            to use only half of that for a batch, so that 2 different jobs (say, AMBER and COBALT)
            can be run in parallel, because using 20 cpus per one job is a waste.
        """
        if ncpus > max_ncpus_per_job:
            # new_ncpus = ncpus
            factor = math.ceil(ncpus / max_ncpus_per_job)
            new_ncpus = ncpus // factor
            # while True:
            #     factor += 1
            #     new_ncpus = ncpus // factor
            #     print(f'ncpus: {ncpus}, factor: {factor}, new_ncpus: {new_ncpus}')
            #     if new_ncpus < max_ncpus_per_job:
            #         print(f'breaking')
            #         break
            if not is_silent:
                logger.info(
                    (msg if msg else 'The number of cpus per batch is ') + f'{ncpus} >{max_ncpus_per_job}. '
                    f'This is usually wasteful, so we are adjusting it '
                    f'to the number <={max_ncpus_per_job}: {new_ncpus} = {ncpus} // {factor}, so '
                    f'{factor} different rules can be run in parallel (say, AMBER and COBALT '
                    f'at the same time).')
            ncpus = new_ncpus
        return ncpus

    if num_batches:
        ncpus_per_batch = max(1, resources.ncpus_available // num_batches)
        resources.ncpus_per_batch = adjust_ncpus_per_job(ncpus_per_batch, max_ncpus_per_job=14, msg=
            f'The number of cpus available is {resources.ncpus_available}, and the number of batches is {num_batches}, '
            f'so the number of cpus per batch would be ')
    if num_samples:
        ncpus_per_sample = max(1, resources.ncpus_available // num_samples)
        resources.ncpus_per_sample = adjust_ncpus_per_job(ncpus_per_sample, max_ncpus_per_job=14, msg=
            f'The number of cpus available is {resources.ncpus_available}, and the number of samples is {num_samples}, '
            f'so the number of cpus per sample would be ')

    if not is_silent:
        if resources.ncpus_local:
            logger.info(f'Local CPUs: {resources.ncpus_local}')
        if ncpus_requested:
            logger.info(f'Total CPUs requested by `umccrise --cores`: {ncpus_requested}')
        logger.info(f'The pipeline can use {resources.ncpus_available} CPUs total.')
        if num_batches:
            logger.info(f'Batches found: {num_batches}, using {resources.ncpus_per_batch} cpus per batch.')
        if num_samples:
            logger.info(f'Samples found: {num_samples}, using {resources.ncpus_per_sample} cpus per sample.')

    return resources


def setup_tibanna(tibanna_id=None, buckets=None):
    try:
        subprocess.check_call(f'tibanna --version', shell=True)
    except subprocess.CalledProcessError:
        logger.err('Error: tibanna is not installed. Please run `pip install -S tibanna`')
        sys.exit(1)

    if not tibanna_id:
        tibanna_id = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) 
                             for _ in range(8))
        assert not check_tibanna_id_exists(tibanna_id), 'Random tibanna ID already exists: ' + tibanna_id

    step_func_name = f'tibanna_unicorn_{tibanna_id}'
    if not check_tibanna_id_exists(tibanna_id):
        buckets_str = '' if not buckets else ('-b ' + ','.join(buckets))
        run_simple(f'tibanna deploy_unicorn -g {step_func_name} {buckets_str} --no-setenv')

    return step_func_name


def check_tibanna_id_exists(tibanna_id):
    step_func_name = f'tibanna_unicorn_{tibanna_id}'
    try:
        subprocess.check_call(f'tibanna list_sfns | grep -w {step_func_name}', shell=True)
    except subprocess.CalledProcessError:
        return False
    else:
        return True
    
    

    