#!/usr/bin/python

'''
Converter from primer sequence files in TXT format into
BED format for amplicon-based mutation calling with Vardict

  example (txt):
chr1    160786646       160786723       ACCAGTTTCTTTCTGAGAACATC GCCTCTTACTGCACTTCATAAAG 1       _CONTROL_V2T4_
chr10   43595874        43595949        CCTTGAAGAAGCCTTATTCTCA  CTGCCTGGTCCACATACAG     4       RET
chr10   43595922        43596043        AGTGGCATTGGGCCTCTA      CGGTACGTGCCGTAGAGAT     1       RET

  example (.primers.bed)
chr1    160786623       160786747       _CONTROL_V2T4_  1       .       160786646       160786724
chr10   43595852        43595969        RET     4       .       43595874        43595950
chr10   43595904        43596063        RET     1       .       43595922        43596044
'''

import sys
import os

if len(sys.argv) != 2:
    print("Usage: " + sys.argv[0] + " <primer_sequence_file.txt>")
    exit(0)

txt_fname = sys.argv[1]
bed_fname = os.path.splitext(txt_fname)[0] + '.primers.bed'

n_lines = 0
with open(txt_fname, 'r') as f:
    with open(bed_fname, 'w') as out:
        for line in f:
            content = line.split('\t')
            chr_name, start, end, primer1_len, primer2_len, x, name = \
                content[0], int(content[1]), int(content[2]), \
                len(content[3]), len(content[4]), int(content[5]), content[6].strip()
            out.write('%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\n' % (chr_name, start - primer1_len, end + primer2_len + 1,
                                                            name, x, '.', start, end + 1))
            n_lines += 1

print("Processed %d lines; Result saved in %s" % (n_lines, bed_fname))
print("Example of conversion:")
with open(txt_fname, 'r') as input:
    with open(bed_fname, 'r') as output:
        print("INPUT (" + txt_fname + "): " + input.readline().strip())
        print("OUTPUT (" + bed_fname + "): " + output.readline().strip())