#!/usr/bin/env python
import os
import sys
from os.path import join, isfile, abspath, dirname, relpath, exists
from setuptools import setup, find_packages

import ngs_utils
from ngs_utils.file_utils import which
from ngs_utils.setup_utils import get_reqs, init, err, compile_tool, find_package_files

name = 'NGS Utils'
package_name = reporting_package_name = ngs_utils.__name__


def install_bedtools():
    success_compilation = compile_tool('BEDtools', ngs_utils.bedtools_dirpath, [ngs_utils.bedtools_execuable_fpath])
    if success_compilation:
        return ngs_utils.bedtools_execuable_fpath
    sys_bedtools_fpath = which('bedtools')
    if sys_bedtools_fpath:
        err('Compilation failed, using bedtools in $PATH: ' + sys_bedtools_fpath)
        return sys_bedtools_fpath
    else:
        err('Compilation of BEDtools at ' + ngs_utils.bedtools_dirpath + ' failed, and no bedtools found in $PATH')
        sys.exit(1)


init(name, package_name, __file__)

print('Installing BEDtools')
bedtools_fpath = install_bedtools()
print('Using BedTools at ' + bedtools_fpath)

setup(
    name=name,
    author='Vlad Saveliev',
    author_email='vlad.saveliev@astrazeneca.com',
    description='Utils for NGS pipelines by Vlad Saveliev',
    keywords='bioinformatics',
    url='https://github.com/vladsaveliev/NGS_Utils',
    license='GPLv3',
    packages=[package_name],
    package_data={
        package_name: [
            'bedtools/*.py',
            'sambamba/*.py',
            relpath(ngs_utils.bedtools_execuable_fpath, abspath(ngs_utils.__name__)),
            relpath(ngs_utils.sambamba_executable_path, abspath(ngs_utils.__name__)),
        ] + find_package_files('reporting', package_name, skip_exts=['.sass', '.coffee'])\
          + find_package_files('reference_data', package_name)
    },
    include_package_data=True,
    zip_safe=False,
    install_requires=get_reqs(),
    setup_requires=['numpy'],
    scripts=[join('scripts', fn) for fn in os.listdir(join(dirname(__file__), 'scripts'))
             if fn.endswith('.py') or fn.endswith('.sh')] +
            [join('scripts', 'converters', fn) for fn in os.listdir(join(dirname(__file__), 'scripts', 'converters'))
             if fn.endswith('.py') or fn.endswith('.sh')],
    classifiers=[
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)

