import os
from datetime import datetime
from os.path import dirname, abspath, join

from ngs_utils import logger
from ngs_utils.file_utils import safe_mkdir


def package_path():
    return dirname(abspath(__file__))


def get_submit_script():
    return join(package_path(), 'submit')


def make_cluster_cmdl(log_dir, app_name=''):
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
        fixed_jobscript = join(log_dir, 'jobscript.sh')
        with open(jobscript) as f_in, open(fixed_jobscript, 'w') as f_out:
            f_out.write(f_in.read().replace('{path}', os.environ["PATH"]))
        cluster_cmdl += f' --jobscript "{fixed_jobscript}"'

    return cluster_cmdl
