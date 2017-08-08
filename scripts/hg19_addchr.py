#!/usr/bin/env python

import sys

inp = sys.stdin
if len(sys.argv) > 1:
    inp = open(sys.argv[1])

out = sys.stdout
if len(sys.argv) > 2:
    out = open(sys.argv[2], 'w')

for l in inp:
    if l and not l.startswith('#') and not l.startswith('chr'):
        l = 'chr' + l
        if l.startswith('chrMT'):
            l = l[:4] + l[5:]
    out.write(l)

inp.close()
out.close()
