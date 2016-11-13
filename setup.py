#!/usr/bin/env python
import os
import sys
from os.path import join, isfile, abspath, dirname
from setuptools import setup, find_packages


from Utils.setup_utils import get_utils_package_files, get_reqs


name = 'Utils'
package_name = reporting_package_name = 'utils'


setup(
    name=name,
    author='Vlad Saveliev',
    author_email='vlad.saveliev@astrazeneca.com',
    description='Vlads bioinf Utils',
    keywords='bioinformatics',
    url='https://github.com/vladsaveliev/Utils',
    license='GPLv3',
    packages=[],
    package_data={
        '.': get_utils_package_files(),
    },
    include_package_data=True,
    zip_safe=False,
    install_requires=get_reqs(),
    setup_requires=['numpy'],
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

