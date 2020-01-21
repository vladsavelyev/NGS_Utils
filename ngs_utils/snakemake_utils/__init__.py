from random import random
import string
import os
import subprocess
import sys
from datetime import datetime
from os.path import dirname, abspath, join
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


def make_cluster_cmdl(log_dir, app_name=''):
    """ Generates cluster command line parameters for snakemake
    """
    from hpc_utils import hpc
    if not hpc.cluster_cmd:
        logger.critical(f'Automatic cluster submission is not supported for the machine "{hpc.name or hpc.hostname}"')

    cluster_submitter = get_submit_script()
    timestamp = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    cluster_cmdl = f' --cluster "{cluster_submitter} {timestamp} {log_dir} {app_name}"'

    # Also overriding jobscript?
    jobscript = hpc.cluster_jobscript
    if jobscript:
        jobscript_file = join(log_dir, 'jobscript.sh')
        with open(jobscript_file, 'w') as f_out:
            f_out.write(jobscript.replace('{path}', os.environ["PATH"]))
        cluster_cmdl += f' --jobscript "{jobscript_file}"'

    return cluster_cmdl


def run_snakemake(snakefile, conf, jobs=None, output_dir=None, forcerun=None,
                  unlock=False, dryrun=False, target_rules=None, cluster=None, cluster_cmd=None,
                  log_dir=None, dag=None, report=None, restart_times=None,
                  tibanna_cfg=None):

    conf['total_cores'] = jobs

    #########################
    #### Setting cluster ####
    #########################

    cluster_param = ''
    cluster_log_dir = ''
    if cluster or cluster_cmd:
        assert log_dir, 'For cluster run, must also specify log_dir'
        if cluster_cmd:
            cluster_param = f' --cluster "{cluster_cmd}"'
        else:
            cluster_log_dir = safe_mkdir(join(log_dir, 'cluster'))
            cluster_param = make_cluster_cmdl(cluster_log_dir, 'umccrise')

    ##########################
    #### Preparing config ####
    ##########################

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
        f'{"--dag " if dag else ""}'
        f'{f"--report {report} " if report else ""}'
        f'{f"--directory {output_dir} " if output_dir else ""}'
        f'{f"-j {jobs} " if jobs else ""}'
        f'--rerun-incomplete '
        f'{f"--restart-times {restart_times - 1}" if restart_times > 1 else ""}'
        f'{cluster_param} '
        f'--configfile {conf_f.name} ' +
        f'{"--dag " if dag else ""}' +
        f'{f"--forcerun {forcerun} " if forcerun else ""}' +
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
        if cluster_log_dir:
            run_simple(f'chmod -R a+r {cluster_log_dir}', silent=True)
            logger.error(f'Review cluster job logs in {cluster_log_dir}')
        sys.exit(1)
    except KeyboardInterrupt:
        logger.error('--------')
        logger.error(f'Interrupted. Fixing logs permissions. Working directory: {output_dir}')
        if cluster_log_dir:
            run_simple(f'chmod -R a+r {cluster_log_dir}', silent=True)
            logger.error(f'Review cluster job logs in {cluster_log_dir}')
        sys.exit(1)
    else:
        logger.info('--------')
        if cluster_log_dir:
            run_simple(f'chmod -R a+r {cluster_log_dir}', silent=True)
        logger.info(f'Finished. Output directory: {output_dir}')


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
    
    
    
    
    
    
    
    
    