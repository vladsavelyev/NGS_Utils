#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

import sys
import os
import shlex
import subprocess
from os.path import basename, join, isfile, getctime

from source.logger import err


def main():
    args = sys.argv[1:]

    if len(args) < 2:
        sys.exit('Usage: ' + __file__ + ' bam sambamba cmdline')

    bam = args[0]
    sambamba = args[1]
    args = args[2:]
    args = [a.replace('__QUOTE__', '"').replace('""', '"') for a in args]
    err(str(args))

    index_bam(bam, sambamba)

    err()
    args = [sambamba] + args
    cmdl = ' '.join((('"' + a + '"') if ' ' in a and not a[0] == '"' else a) for a in args)
    err(cmdl)
    ret_code = subprocess.call(cmdl, shell=True)
    if ret_code != 0:
        err()
        err('Ret code = ' + str(ret_code) + ', retrying...')
        indexed_bam = bam + '.bai'
        if isfile(indexed_bam):
            os.remove(indexed_bam)
        index_bam(bam, sambamba)
        subprocess.call(cmdl, shell=True)


def index_bam(bam_fpath, sambamba):
    indexed_bam = bam_fpath + '.bai'
    if isfile(indexed_bam):
        return
    # if not isfile(indexed_bam) or getctime(indexed_bam) < getctime(bam_fpath):
    err('Indexing BAM, writing ' + indexed_bam + '...')
    cmdline = '{sambamba} index {bam_fpath}'.format(**locals())
    subprocess.call(cmdline, shell=True)


if __name__ == '__main__':
    main()
