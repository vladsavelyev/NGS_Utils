#!/usr/bin/python

'''
Converter from BED + primer sequence files in TXT format into
traditional 8-column BED format

  example (BEDfile.txt, 1-based coordinates):
chr13   32890471        32890723        BRCA2_ex02_01
chr17   41276004        41276264        BRCA1_ex02_01
  example (Primerfile.txt):
BRCA2_ex02_01-F AAGACTCGGCAGCATCTCCACAGGAGATGGGACTGAATTAGAATTCAAAC
BRCA1_ex02_01-F AAGACTCGGCAGCATCTCCATTGCATAGGAGATAATCATAGGAATC
BRCA2_ex02_01-R GCGATCGTCACTGTTCTCCAACACTGTGACGTACTGGGTTTTTAGCAAGC
BRCA1_ex02_01-R GCGATCGTCACTGTTCTCCAGAATCTTTAAAAATAAAGGACGTTGTCATT

  example (output.primers.bed, 0-based coordinates)
chr13   32890420    32890773    BRCA2_ex02_01   .   .   32890470    32890723
chr17   41275957    41276314    BRCA1_ex02_01   .   .   41276003    41276264
'''

import sys
import os

if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + " <Regions File.txt> <Primer Sequence File.txt> [output.bed]")
    exit(0)

txt_bed_fname = sys.argv[1]
txt_primer_fname = sys.argv[2]
if len(sys.argv) > 3:
    out_bed_fname = sys.argv[3]
else:
    out_bed_fname = os.path.splitext(txt_bed_fname)[0] + '.primers.bed'

primers = dict()
with open(txt_primer_fname, 'r') as f:
    for line in f:
        primers[line.split()[0]] = len(line.split()[1])

n_lines = 0
with open(txt_bed_fname, 'r') as f:
    with open(out_bed_fname, 'w') as out:
        for line in f:
            content = line.split()
            if not content:
                break
            chr_name, start, end, name = content[0], int(content[1]), int(content[2]), content[3].strip()
            if name + '-F' not in primers.keys() or name + '-R' not in primers.keys():
                sys.stderr.write('Error: no corresponding primer sequences for %s\n' % name)
                sys.exit(1)
            primer1_len = primers[name + '-F']
            primer2_len = primers[name + '-R']
            out.write('%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\n' % (chr_name, start - primer1_len - 1, end + primer2_len,
                                                            name, '.', '.', start - 1, end))
            n_lines += 1

print("Processed %d lines and %d primer sequences (F and R); Result saved in %s" % (n_lines, len(primers), out_bed_fname))
print("Example of conversion:")
with open(txt_bed_fname, 'r') as input:
    with open(out_bed_fname, 'r') as output:
        print("INPUT (" + txt_bed_fname + "): " + input.readline().strip())
        print("OUTPUT (" + out_bed_fname + "): " + output.readline().strip())
