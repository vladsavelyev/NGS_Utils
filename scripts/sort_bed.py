#!/usr/bin/env python
from os.path import basename

import sys
from optparse import OptionParser

from ngs_utils.file_utils import adjust_path
from ngs_utils.logger import critical
from ngs_utils.bed_utils import verify_bed, sort_bed


def main():
    parser = OptionParser(usage='Usage: ' + basename(__file__) + ' -o Output_BED_file -g hg19 Input_BED_file')
    parser.add_option('-o', '--output-bed', dest='output_fpath')
    parser.add_option('-g', '--genome', dest='genome')
    parser.add_option('--fai', dest='fai_fpath')
    (opts, args) = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    if not opts.output_fpath:
        critical(parser.usage)

    sort_bed(input_bed_fpath=verify_bed(args[0], is_critical=True),
             output_bed_fpath=adjust_path(opts.output_fpath),
             fai_fpath=adjust_path(opts.fai_fpath),
             genome=opts.genome)


if __name__ == '__main__':
    main()
