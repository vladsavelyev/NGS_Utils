#!/usr/bin/env python
import sys
py_v = sys.version_info[:2]
if not (py_v == (2, 7) or py_v >= (3, 3)):
    sys.exit('Only Python 2.7 or 3.3 and up are supported. Current version: ' + '.'.join(py_v))
    
from os.path import join
from ngs_utils import setup_utils


import ngs_utils
package_name = ngs_utils.__name__


version = setup_utils.init(package_name, package_name, __file__)


from setuptools import setup
setup(
    name=package_name,
    version=version,
    author='Vlad Saveliev',
    author_email='vlad.saveliev@astrazeneca.com',
    description='Utils for NGS pipelines by Vlad Saveliev',
    keywords='bioinformatics',
    url='https://github.com/vladsaveliev/NGS_Utils',
    license='GPLv3',
    packages=[
        package_name,
        'tab_utils'
    ],
    package_data={
        package_name: setup_utils.find_package_files('', package_name, skip_exts=['.sass', '.coffee'])
    },
    include_package_data=True,
    zip_safe=False,
    install_requires=setup_utils.get_reqs(),
    setup_requires=['numpy'],
    scripts=[path for path in
         [join('scripts', fn) for fn in [
             'tabutils',
             'standardize_bed.py',
             'split_bed.py',
             'qstat.py',
             'sort_bed.py',
             'bed_file_from_gene_list.py',
             'html_to_base64.py',
             'group_and_merge_by_gene.py'
         ]]],
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
