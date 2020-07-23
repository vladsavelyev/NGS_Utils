#!/usr/bin/env python
from os.path import join
from setuptools import setup

import ngs_utils
package_name = ngs_utils.__name__

package_data = None
try:
    import versionpy
except ImportError:
    version = 'dev'
else:
    version = versionpy.get_version(package_name)
    package_data = {
        package_name: versionpy.find_package_files('', package_name, skip_exts=['.sass', '.coffee'])
    }
print('package_data is None:', package_data is None)


setup(
    name=package_name,
    version=version,
    author='Vlad Saveliev',
    author_email='vladislav.sav@gmail.com',
    description='Utils for NGS pipelines by Vlad Saveliev',
    keywords='bioinformatics',
    url='https://github.com/vladsaveliev/NGS_Utils',
    license='GPLv3',
    packages=[package_name],
    package_data=package_data,
    include_package_data=True,
    zip_safe=False,
    install_requires=['click'],
    setup_requires=['numpy'],
    scripts=[path for path in
         [join('scripts', fn) for fn in [
             'standardize_bed.py',
             'split_bed.py',
             'qstat.py',
             'sort_bed.py',
             'bed_file_from_gene_list.py',
             'html_to_base64.py',
             'group_and_merge_by_gene.py',
             'hg19_addchr.py',
             'generate_bed.py',
         ]]],
    classifiers=[
        'Environment :: Console',
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
