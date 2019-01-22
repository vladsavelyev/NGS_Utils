#!/usr/bin/env python3.3

import gzip, sys

if __name__=='__main__':
    feature = sys.argv[1] # has to be either 'cds' or 'exon' or 'gene'

    gene_list = None
    feature_in_gene_list = set([])
    if len(sys.argv) == 3:
        gene_list = set([])
        with open(sys.argv[2]) as gl:
            for line in gl:
                toks = line.rstrip().split('\t')
                gene_list.add( toks[0] )

    for ln, line in enumerate(gzip.open('refGene.txt.gz', 'rb')):
        toks = line.decode().rstrip().split('\t')
        _, tid, chrom, strand, tx_start, tx_end, cds_start, cds_end, exon_count, exon_starts, exon_ends, _, symbol, _, _, _ = toks
    
        chrom = chrom.replace('chr','')
        if chrom.find('_') >= 0:
            continue
    
        tx_start, tx_end = int(tx_start), int(tx_end)
        cds_start, cds_end = int(cds_start), int(cds_end)
        exon_count = int(exon_count)
        exon_starts = [ int(v) for v in exon_starts.split(',') if v ]
        exon_ends = [ int(v) for v in exon_ends.split(',') if v ]

        if gene_list is not None:
            if symbol not in gene_list:
                continue
            else:
                feature_in_gene_list.add( symbol )

   
        if feature == 'exon': 
            for s, e in zip(exon_starts, exon_ends):
                print( '{}\t{}\t{}\t{}|{}'.format(chrom, s, e, symbol, tid) )
        elif feature == 'cds':
            # non-coding
            if cds_start == cds_end:
                continue

            for s, e in zip(exon_starts, exon_ends):
                if s <= cds_start <= e:
                    s = cds_start
                if s <= cds_end <= e:
                    e = cds_end
                print('{}\t{}\t{}\t{}|{}'.format(chrom, s, e, symbol, tid))
        elif feature == 'gene':
            print('{}\t{}\t{}\t{}'.format(chrom, min(exon_starts), max(exon_ends), symbol) )

  
    if gene_list is not None:
        print('Checking for genes NOT INCLUDED in refGene', file=sys.stderr) 
        for g in sorted(list(gene_list)) :
            if g not in feature_in_gene_list:
                print('Gene in list that is NOT in refGene: {}'.format( g ), file=sys.stderr )
