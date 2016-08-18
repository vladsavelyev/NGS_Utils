from os.path import dirname, join, abspath, isdir, isfile
import os
from sys import platform as sys_platform
import platform


bedtools_dirpath = abspath(join(dirname(__file__), 'bedtools', 'bedtools2'))
bedtools_execuable_fpath = join(bedtools_dirpath, 'bin', 'bedtools')


sambamba_bin_dirpath = abspath(join(dirname(__file__), 'sambamba', 'bin'))
if 'darwin' in sys_platform:
    sambamba_executable_path = join(sambamba_bin_dirpath, 'sambamba_osx')
elif 'redhat' in platform.dist() or 'centos' in platform.dist():
    sambamba_executable_path = join(sambamba_bin_dirpath, 'sambamba_centos')
else:
    sambamba_executable_path = join(sambamba_bin_dirpath, 'sambamba_lnx')