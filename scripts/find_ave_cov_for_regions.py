#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

import sys
from logging import info

from source.file_utils import adjust_path
from source.targetcov.Region import save_regions_to_bed, GeneInfo


def main():
    args = sys.argv[1:]

    if len(args) < 2:
        sys.exit('Usage: ' + __file__ + ' sambamba_depth_report sample_name bed_col_num')

    bedcov_hist_fpath, sample_name, bed_col_num = args

    amplicons = summarize_bedcoverage_hist_stats(adjust_path(bedcov_hist_fpath), sample_name, int(bed_col_num))

    amplicons = sorted(amplicons, key=lambda a: (a.chrom, a.gene_name, a.start))

    for r in amplicons:
        r.calc_avg_depth()

    save_regions_to_seq2cov_output__nocnf(sample_name, amplicons)


def save_regions_to_seq2cov_output__nocnf(sample_name, regions, output_fpath=None):
    final_regions = []
    gene = None
    total_amplicons = 0
    total_genes = 0
    for a in regions:
        a.sample_name = sample_name
        a.feature = 'Amplicon'
        if gene and gene.gene_name != a.gene_name:
            __sum_up_gene(gene)
            final_regions.append(gene)
            total_genes += 1
            gene = None
        if not gene:
            gene = GeneInfo(sample_name=sample_name, gene_name=a.gene_name, chrom=a.chrom,
                            strand=a.strand, feature='Whole-Gene')
        gene.add_amplicon(a)
        final_regions.append(a)
        total_amplicons += 1
    if gene:
        __sum_up_gene(gene)
        final_regions.append(gene)

    info('Number of final regions: ' + str(len(final_regions)) + ', out of them - ' +
         str(total_amplicons) + ' amplicons, ' + str(total_genes) + ' genes.')
    coverage_info = []
    for r in final_regions:
       # if r.avg_depth is not None:  #and r.avg_depth != 0:
        coverage_info.append([sample_name, r.gene_name, r.chrom, r.start + 1, r.end, r.feature, r.size, r.avg_depth])

    info('Coverage info lines: ' + str(len(coverage_info)))
    if output_fpath:
        with open(output_fpath, 'w') as f:
            for fs in coverage_info:
                f.write('\t'.join(map(str, fs)) + '\n')
    else:
        for fs in coverage_info:
            print '\t'.join(map(str, fs))


def __sum_up_gene(g):
    g.start = g.amplicons[0].start
    g.end = g.amplicons[-1].end
    g.size = sum(a.end - a.start for a in g.amplicons)
    g.avg_depth = sum(float(a.size) * a.avg_depth for a in g.amplicons) / g.get_size() if g.get_size() else 0


if __name__ == '__main__':
    main()