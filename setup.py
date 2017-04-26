#!/usr/bin/env python
import os
import sys
from os.path import join, isfile, abspath, dirname, relpath, exists


os.subprocess.call('pip install --upgrade pip', shell=True)
os.subprocess.call('pip install --upgrade --ignore-installed setuptools', shell=True)


import ngs_utils
from ngs_utils.setup_utils import get_reqs, init, err, find_package_files

package_name = ngs_utils.__name__
print package_name


init(package_name, package_name, __file__)
print('')


from setuptools import setup
setup(
    name=package_name,
    version='1.0',
    author='Vlad Saveliev',
    author_email='vlad.saveliev@astrazeneca.com',
    description='Utils for NGS pipelines by Vlad Saveliev',
    keywords='bioinformatics',
    url='https://github.com/vladsaveliev/NGS_Utils',
    license='GPLv3',
    packages=[package_name, 'tab_utils'],
    package_data={
        package_name: [
        ] + find_package_files('reporting', package_name, skip_exts=['.sass', '.coffee'])
          + find_package_files('reference_data', package_name)
    },
    include_package_data=True,
    zip_safe=False,
    install_requires=get_reqs(),
    setup_requires=['numpy'],
    scripts=[path for path in
       [join('scripts', fn) for fn in os.listdir(join(dirname(__file__), 'scripts'))] +
       [join('scripts', 'converters', fn) for fn in os.listdir(join(dirname(__file__), 'scripts', 'converters'))]
       if isfile(path)],
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
