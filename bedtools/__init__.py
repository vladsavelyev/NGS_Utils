from os.path import dirname, join, abspath, isdir, isfile
import os

from Utils.file_utils import file_exists, which
from Utils.logger import critical, err
from Utils import bedtools_execuable_fpath


def find_executable():
    if not file_exists(bedtools_execuable_fpath):
        exec_fpath_in_path = which('bedtools')
        if exec_fpath_in_path:
            err('BEDtools compilation failed, using executable found in $PATH: ' + exec_fpath_in_path + '\n')
        else:
            critical('Error: could not find BEDtools executable at ' + bedtools_execuable_fpath
                     + ' or in $PATH')
    return bedtools_execuable_fpath


if not isfile(bedtools_execuable_fpath):
    critical('Error: BEDtools execuatable not found under ' + bedtools_execuable_fpath +
             '! Please first install BEDtools, then import this module.')
os.environ['PATH'] = dirname(bedtools_execuable_fpath) + ':' + os.environ['PATH']
# noinspection PyUnresolvedReferences
from pybedtools import BedTool
