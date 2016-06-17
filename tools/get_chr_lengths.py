#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

import sys

from source.logger import err, critical
from source.file_utils import file_exists, verify_file, file_transaction
from source.utils import get_chr_lengths_from_seq


def main():
    if len(sys.argv) <= 2:
        critical('Usage: ' + __file__ + ' path_to_.fa')

    seq_fpath = sys.argv[1]
    seq_fpath = verify_file(seq_fpath, is_critical=True)
    chr_lengths = get_chr_lengths_from_seq(seq_fpath)

    for c, l in chr_lengths:
        sys.stdout.write(c + '\t' + str(l) + '\n')


if __name__ == '__main__':
    main()
