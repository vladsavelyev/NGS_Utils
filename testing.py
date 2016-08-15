import subprocess
import unittest
import os
from genericpath import getsize, getmtime
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath

import shutil
import sys


def info(msg=''):
    sys.stdout.write(msg + '\n')

def call(cmdl):
    info(cmdl)
    subprocess.call(cmdl)

def check_call(cmdl):
    info(cmdl if isinstance(cmdl, basestring) else ' '.join(cmdl))
    subprocess.check_call(cmdl, shell=isinstance(cmdl, basestring))


class BaseTestCase(unittest.TestCase):
    script = None

    data_dir = 'results'
    results_dir = 'results'
    gold_standard_dir = 'gold_standard'

    remove_work_dir_on_success = False

    def setUp(self):
        if not isdir(self.data_dir):
            os.makedirs(self.data_dir)
        if not exists(self.results_dir):
            os.makedirs(self.results_dir)

    def _check_file(self, fpath, diff_ignore_re=''):
        assert isfile(fpath)
        assert getsize(fpath) > 0
        if isdir(self.gold_standard_dir):
            cmp_fpath = join(self.gold_standard_dir, relpath(fpath, self.results_dir))
            if cmp_fpath and isfile(cmp_fpath):
                if diff_ignore_re:
                    cmdl = ['diff', '-q', '--ignore-matching-lines', diff_ignore_re, fpath, cmp_fpath]
                else:
                    cmdl = ['diff', '-q', fpath, cmp_fpath]
                check_call(cmdl)
