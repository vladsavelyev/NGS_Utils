#!/usr/bin/env python
import os, glob, sys, shutil, subprocess
from os.path import join, isfile, abspath, dirname, relpath, isdir
import click
import importlib


VERSION_COMPONENTS = ['MAJOR', 'MINOR', 'BUGFIX']

def validate_version(ctx, param, value):
    if '.' in value:
        if len(value.split('.')) != 3:
            raise click.BadParameter(f'Version must have 3 components. Got: {value}')
    else:
        value = value.upper()
        if value not in VERSION_COMPONENTS:
            raise click.BadParameter(f'Parameter must be either a 3-component version tag, or one of {VERSION_COMPONENTS}')
    return value


@click.command('release')
@click.argument('new_version', default='BUGFIX', callback=validate_version)
@click.option('-p', 'package_name')
def release(new_version, package_name=None):
    ''' Usage: release [bugfix,minor,major,1.0.1] [-p package_name]
    '''
    version_file, new_version = increment_version(new_version, package_name)
    run_cmdl(f'git add {version_file}')
    run_cmdl(f'git commit -m "Bump {new_version}"')
    run_cmdl(f'git tag {new_version}')
    run_cmdl(f'git push')
    run_cmdl(f'git push --tags')


def err(msg=''):
    sys.stderr.write(msg + '\n')


def get_git_revision():
    try:
        import subprocess
        git_revision = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).rstrip()
    except:
        git_revision = ''
        pass
    if isinstance(git_revision, bytes):
        git_revision = git_revision.decode()
    return git_revision


def get_cur_version(package_name='*'):
    version_py = glob.glob(f'{package_name}/_version.py')
    if version_py:
        version_py = version_py[0]
        package_name = dirname(version_py)
        cur_version = importlib.import_module(f'{package_name}._version').__version__
        err(f'Found current version file {version_py}, current version: {cur_version}')
    else:
        version_txt = 'VERSION.txt'
        if isfile(version_txt):
            cur_version = open(version_txt).read().strip()
            err(f'Could not find {package_name}/_version.py, but found {version_txt} - read version {cur_version}')
        else:
            cur_version = '0.1.0'
            err(f'Could not find {package_name}/_version.py, inititalising with version {cur_version}')
    return cur_version


def increment_version(arg='BUGFIX', package_name=None):
    package_name = package_name or '*'

    if arg in VERSION_COMPONENTS:
        cur_version = get_cur_version(package_name)
        components = cur_version.split('.')
        assert len(components) == 3, f'Version must have 3 components. Cannot parse "{cur_version}"'

        comp_ind = VERSION_COMPONENTS.index(arg)
        err(f'Incrementing {arg} component: "{components[comp_ind]}"')
        components[comp_ind] = str(int(components[comp_ind]) + 1)
        new_version = '.'.join(components)

    else:
        new_version = arg

    version_py = glob.glob(f'{package_name}/_version.py')
    if version_py:
        version_py = version_py[0]
    else:
        if package_name == '*':
            err(f'Please specify package_name to initiate version file package_name/_version.py')
            sys.exit(1)
        else:
            version_py = f'{package_name}/_version.py'
            err(f'Could not find current version file under {version_py}, creating a new file {version_py}')

    git_rev = get_git_revision()
    with open(version_py, 'w') as f:
        f.write((
            '# Do not edit this file, pipeline versioning is governed by git tags\n' +
            '__version__ = \'' + new_version + '\'\n' +
            '__git_revision__ = \'' + str(git_rev) + '\'') + '\n')

    err(f'New version: {new_version}, written to {version_py}')
    return version_py, new_version


# def init(package_name, setup_py_fpath):
#     if abspath(dirname(setup_py_fpath)) != abspath(os.getcwd()):
#         sys.exit('Please, change to ' + dirname(setup_py_fpath) + ' before running setup.py\n')
#     return get_cur_version(package_name)


def clean_package(package_name, dirpath='.'):
    print('Cleaning up binary, build and dist for ' + package_name + ' in ' + dirpath + '...')
    if isdir(join(dirpath, 'build')):
        shutil.rmtree(join(dirpath, 'build'))
    if isdir(join(dirpath, 'dist')):
        shutil.rmtree(join(dirpath, 'dist'))
    if isdir(join(dirpath, package_name + '.egg-info')):
        shutil.rmtree(join(dirpath, package_name + '.egg-info'))
    print('Done.')


def get_reqs():
    try: # for pip >= 10
        from pip._internal.req import parse_requirements
    except ImportError: # for pip <= 9.0.3
        from pip.req import parse_requirements
    
    try:
        install_reqs = parse_requirements('requirements.txt', session=False)
    except TypeError:
        install_reqs = parse_requirements('requirements.txt')
    reqs = [str(ir.req) for ir in install_reqs if ir.req]
    return reqs


def find_package_files(dirpath, package, skip_exts=None):
    paths = []
    for (path, dirs, fnames) in os.walk(join(package, dirpath)):
        for fname in fnames:
            if skip_exts and any(fname.endswith(ext) for ext in skip_exts):
                continue
            fpath = join(path, fname)
            paths.append(relpath(fpath, package))
    return paths


# ''' Versioning:
# 1. Write each version to VERSION.txt
# 2. If the changes are significant, tag the release and push the new tag:
#    $ python setup.py tag '''
# def write_version_py(package_name):
#     version_txt = glob.glob('*/VERSION.txt') or glob.glob('VERSION.txt')
#     if not version_txt:
#         return None
#
#     with open(version_txt) as f:
#         v = f.read().strip().split('\n')[0]
#
#     try:
#         import subprocess
#         git_revision = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).rstrip()
#     except:
#         git_revision = ''
#         pass
#     if isinstance(git_revision, bytes):
#         git_revision = git_revision.decode()
#
#     version_py = os.path.join(package_name, 'version.py')
#     with open(version_py, 'w') as f:
#         f.write((
#             '# Do not edit this file, pipeline versioning is governed by git tags\n' +
#             '__version__ = \'' + v + '\'\n' +
#             '__git_revision__ = \'' + str(git_revision) + '\'') + '\n')
#     return v


def run_cmdl(_cmd):
    print('$ ' + _cmd)
    subprocess.run(_cmd, shell=True, check=True)


def compile_tool(tool_name, dirpath, requirements):
    if not all(isfile(join(dirpath, req)) for req in requirements):
        print('Compiling ' + tool_name)
        run_cmdl('make -C ' + dirpath)
        if not all(isfile(join(dirpath, req)) for req in requirements):
            err('Failed to compile ' + tool_name + ' (' + dirpath + ')\n')
            return False
    return True


def is_cleaning():
    return len(sys.argv) == 2 and sys.argv[1] == 'clean'
