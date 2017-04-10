#!/usr/bin/env python
from __future__ import print_function
from __future__ import print_function
from __future__ import print_function
import sys
from optparse import OptionParser
import subprocess

from ngs_utils.bed_venn.venn import run, save_venn_diagram_data, write_html
from os.path import dirname, basename, join, splitext, isfile
import json
import itertools

from ngs_utils.bed_utils import verify_bed
from ngs_utils.file_utils import adjust_path, safe_mkdir
from ngs_utils.logger import critical


beds_dirpath = '/Users/vlad/googledrive/az/reference_data/genomes/Hsapiens/hg19/bed/intersection/merged'


def main():
    parser = OptionParser(usage='Usage: ' + basename(__file__) + ' bed1 bed2 ... -o results_dir')
    parser.add_option('-o', '--output-dir', dest='output_dir')
    (opts, args) = parser.parse_args(sys.argv[1:])

    if len(args) == 0:
        bed_fpaths = tuple([join(beds_dirpath, bed) for bed in [
            'AZ.bed',
            'CRE.bed',
            'IDT.bed',
            'Med.bed',
            # 'One.bed',
            # 'One.hg38tohg19.bed',
            'V6.bed',
        ]])
        names_map = {
            'AZ': 'AZ Exome',
            'CRE': 'B',
            'V6': 'A',
            'IDT': 'D',
            'Med': 'C'}
    else:
        bed_fpaths = [verify_bed(bed) for bed in args]
        names_map = dict()

    if not opts.__dict__.get('output_dir'):
        output_dir = beds_dirpath
    else:
        output_dir = adjust_path(opts.output_dir)
        safe_mkdir(output_dir)
    work_dirpath = safe_mkdir(join(output_dir, 'intersections'))

    intersection_size_by_subset = run(work_dirpath, bed_fpaths)

    json_txt = save_venn_diagram_data(intersection_size_by_subset, names_map)
    
    html_file = write_html(output_dir, json_txt, bed_fpaths)

    print('-----------------------')
    print('')
    print('HTML: ' + html_file)


if __name__ == '__main__':
    main()


