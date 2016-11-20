import subprocess
import unittest
import os
from genericpath import getsize, getmtime
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath, realpath
import shutil
import sys
from datetime import datetime

from ngs_utils.file_utils import verify_dir, verify_file


def info(msg=''):
    sys.stdout.write(msg + '\n')

def call(cmdl):
    info(cmdl)
    subprocess.call(cmdl)

def check_call(cmdl):
    info(cmdl if isinstance(cmdl, basestring) else ' '.join(cmdl))
    subprocess.check_call(cmdl, shell=isinstance(cmdl, basestring))

def swap_output(output_path):
    if not exists(output_path):
        return

    last_changed = datetime.fromtimestamp(getmtime(output_path))
    prev_output_path = output_path + '_' + last_changed.strftime('%Y_%m_%d__%H_%M_%S')
    os.rename(output_path, prev_output_path)
    return prev_output_path

def swap_prev_symlink(output_path, prev_output_path):
    prev_output_link = get_prev(output_path)
    if exists(prev_output_link):
        os.remove(prev_output_link)
    os.symlink(prev_output_path, prev_output_link)

def get_prev(fpath):
    return fpath + '_prev'


class BaseTestCase(unittest.TestCase):
    script = None

    data_dir = 'data'
    results_dir = 'results'
    gold_standard_dir = 'gold_standard'

    remove_work_dir_on_success = False

    def setUp(self):
        if not isdir(self.data_dir):
            os.makedirs(self.data_dir)
        if not exists(self.results_dir):
            os.makedirs(self.results_dir)

    def _check_file(self, fpath, diff_ignore_re=''):
        assert isfile(fpath), fpath
        assert getsize(fpath) > 0, fpath

        cmp_fpath = None
        if isdir(self.gold_standard_dir):
            cmp_fpath = join(self.gold_standard_dir, relpath(fpath, self.results_dir))
        elif isfile(get_prev(fpath)):
            cmp_fpath = get_prev(fpath)

        if cmp_fpath and isfile(cmp_fpath):
            if diff_ignore_re:
                cmdl = ['diff', '-q', '--ignore-matching-lines', diff_ignore_re, fpath, cmp_fpath]
            else:
                cmdl = ['diff', '-q', fpath, cmp_fpath]
            check_call(cmdl)


    def _check_dir_not_empty(self, dirpath, description=None):
        assert verify_dir(dirpath, description=description), dirpath
        contents = [join(dirpath, fname) for fname in os.listdir(dirpath)
                    if not fname.startswith('.')]
        assert len(contents) >= 1, dirpath + ': ' + str(contents)
        assert all(verify_file(realpath(fpath), is_critical=True)
                   for fpath in contents
                   if isfile(realpath(fpath))), dirpath + ': ' + str(contents)
