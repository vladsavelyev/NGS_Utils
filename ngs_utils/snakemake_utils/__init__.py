import os
from datetime import datetime
from os.path import dirname, abspath, join
from ngs_utils import logger
from ngs_utils.file_utils import safe_mkdir
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
    from hpc_utils.hpc import get_loc
    loc = get_loc()
    if not loc.cluster:
        logger.critical(f'Automatic cluster submission is not supported for the machine "{loc.name}"')

    cluster_submitter = get_submit_script()
    timestamp = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    cluster_cmdl = f' --cluster "{cluster_submitter} {timestamp} {log_dir} {app_name}"'

    # Also overriding jobscript?
    jobscript = loc.cluster.get('jobscript')
    if jobscript:
        jobscript_file = join(log_dir, 'jobscript.sh')
        with open(jobscript_file, 'w') as f_out:
            f_out.write(jobscript.replace('{path}', os.environ["PATH"]))
        cluster_cmdl += f' --jobscript "{jobscript_file}"'

    return cluster_cmdl


def run_snakemake(smk_file, conf, jobs=None, output_dir=None, force_rerun=None,
                  unlock=False, dryrun=False, target_rules=None):
    """ Runs snakemake
    """
    f = tempfile.NamedTemporaryFile(mode='wt', delete=False)
    yaml.dump(conf, f)
    f.close()

    if isinstance(target_rules, str):
        target_rules = [target_rules]

    cmd = (f'snakemake ' +
           f'{" ".join(target_rules) if target_rules else ""} ' +
           f'--snakefile {smk_file} ' +
           f'--printshellcmds ' +
           f'{"--dryrun " if dryrun else ""}' +
          (f'--directory {output_dir} ' if output_dir else '') +
          (f'-j {jobs} ' if jobs is not None else '') +
           f'--rerun-incomplete ' +
           f'--configfile {f.name} ' +
          (f'--forcerun {force_rerun}' if force_rerun else '')
           )

    if unlock:
        print('* Unlocking previous run... *')
        run_simple(cmd + ' --unlock')
        print('* Now rerunning *')
    run_simple(cmd)

