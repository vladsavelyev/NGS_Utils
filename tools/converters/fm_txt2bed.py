#!/usr/bin/python

'''
Converter from Foundation Medicine TXT format into BED format

  example: input (txt):
Target  interval
WG_snp_1        chr1:852745-852884
TNFRSF14_target_1       chr1:2488073-2488202
TNFRSF14_target_2       chr1:2489154-2489283

  example: output (bed)
chr1    852745       852884       WG_snp_1
chr1    2488073      2488202      TNFRSF14_target_1
chr1    2489154      2489283      TNFRSF14_target_2
'''

import sys
import os

if len(sys.argv) != 2:
    print("Usage: " + sys.argv[0] + " <input.txt>")
    exit(0)

txt_fname = sys.argv[1]
bed_fname = os.path.splitext(txt_fname)[0] + '.bed'

n_lines = 0
example_line = ''
with open(txt_fname, 'r') as f:
    with open(bed_fname, 'w') as out:
        for line in f:
            if ':' not in line and '-' not in line:  # comment line
                continue
            elif not example_line:
                example_line = line
            target, interval = line.split('\t')[0], line.split('\t')[1]
            chr_name = interval.split(':')[0]
            start, end = int(interval.split(':')[1].split('-')[0]), int(interval.split(':')[1].split('-')[1])
            out.write('%s\t%d\t%d\t%s\n' % (chr_name, start, end, target))
            n_lines += 1

print("Processed %d lines; Result saved in %s" % (n_lines, bed_fname))
print("Example of conversion:")
with open(bed_fname, 'r') as output:
    print("INPUT (" + txt_fname + "): " + example_line.strip())
    print("OUTPUT (" + bed_fname + "): " + output.readline().strip())
