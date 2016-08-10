from os.path import dirname, join, abspath, isdir, isfile
import os

from Utils.file_utils import file_exists
from Utils.logger import critical


def find_executable():
    exec_fpath = abspath(join(dirname(__file__), 'bedtools2', 'bin', 'bedtools'))
    if not file_exists(exec_fpath):
        critical('Error: could not find BedTools executable at ' + exec_fpath)
    return exec_fpath


os.environ['PATH'] = dirname(find_executable()) + ':' + os.environ['PATH']
# noinspection PyUnresolvedReferences
from pybedtools import BedTool
