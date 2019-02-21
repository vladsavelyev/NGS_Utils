#!/usr/bin/env python
import os
import sys
from bed_annotation import canon_transcript_per_gene
from ngs_utils.file_utils import open_gzipsafe
from ngs_utils.key_genes_utils import get_genes_from_file
from ngs_utils.logger import warn, critical
from ngs_utils.reference_data import get_key_genes_set
import click

f''' * Examples *
     
     Generate BED files for key genes:

        python {__file__} -g GRCh37 --key-genes --features CDS,stop_codon\
            > umccr_cancer_genes.GRCh37.coding.bed
        python {__file__} -g GRCh37 --key-genes --features transcript\
            > umccr_cancer_genes.GRCh37.transcript.bed

     Generate coding regions for SAGE (https://github.com/hartwigmedical/hmftools/tree/master/sage)

        python {__file__} -g GRCh37 --principal --key-genes --biotypes protein_coding,decay --features CDS\
           | sort -k1,1V -k2,2n | grep -v ^MT | grep -v ^GL\
           | bedtools merge -c 4 -o collapse -i -\
           > coding_regions.bed
'''

@click.command()
@click.option('-g', 'genome', default='GRCh37')
@click.option('--gtf', 'gtf_path', help='Path to GTF to extract features. Default is by hpc_utils')
@click.option('--all-transcripts', 'all_transcripts', is_flag=True, help='Use all transcripts (default is principal+alternative by APPRIS)')
@click.option('-p', '--principal', 'principal', is_flag=True, help='Use only principal transcripts (by APPRIS)')
@click.option('-k', '--key-genes', 'only_key_genes', is_flag=True, help='Use UMCCR key cancer genes only')
@click.option('--genes', 'gene_list', help='Use genes from the list only')
@click.option('--biotypes', 'biotypes', default='protein_coding,decay', help='Feature types to extract')
@click.option('--features', 'features', default='CDS,stop_codon', help='Feature types to extract')

def main(genome=None, gtf_path=None, all_transcripts=False, principal=False, only_key_genes=False, gene_list=None,
         biotypes='', features=''):
    out = sys.stdout

    # GTF
    if not gtf_path:
        try:
            from hpc_utils import hpc
        except ImportError:
            critical('GTF file is needed. Either install hpc_utils, or provide GTF with --gtf')
        else:
            if genome == 'GRCh37':
                gtf_path = os.path.join(hpc.get_ref_file(key='pyensembl_data'), 'GRCh37/ensembl75/Homo_sapiens.GRCh37.75.gtf.gz')
            else:
                gtf_path = os.path.join(hpc.get_ref_file(key='pyensembl_data'), 'GRCh38/ensembl95/Homo_sapiens.GRCh38.95.gtf.gz')

    # Genes
    key_genes = None
    if only_key_genes:
        key_genes = get_key_genes_set()
    elif gene_list:
        key_genes = get_genes_from_file(gene_list)

    # Transcripts
    transcripts_by_gid = None
    if not all_transcripts:
        if principal:
            transcripts_by_gid = {
                k: [v] for k, v in
                canon_transcript_per_gene(genome, only_principal=True, use_gene_id=True).items()
            }
        else:
            transcripts_by_gid = canon_transcript_per_gene(genome, use_gene_id=True)

    # Options
    biotypes = biotypes.strip()
    if biotypes:
        biotypes = biotypes.split(',')
    features = features.strip()
    if features:
        features = features.split(',')

    genes_set = set()
    genes_without_canon = set()
    warn(f'Parsing {gtf_path}')
    with open_gzipsafe(gtf_path) as f:
        lines_cnt = 0
        region_cnt = 0
        for l in f:
            if not l.startswith('#') and l.strip():
                lines_cnt += 1
                fields = l.strip().split('\t')
                try:
                    chrom, biotype, feature, start, end, _, _, _, annotations = fields
                except:
                    warn(f'Cannot read fields {str(fields)}')
                    raise

                if biotypes:
                    if not any(bt in biotype for bt in biotypes):
                        continue

                if features:
                    if not any(feature == ft for ft in features):
                        continue

                annotations = {kv.split()[0].strip().strip('"'):
                                   kv.split()[1].strip().strip('"')
                               for kv in annotations.split('; ')}

                gene_name = annotations['gene_name']
                if only_key_genes and gene_name not in key_genes:
                    continue

                if transcripts_by_gid:
                    gene_id = annotations['gene_id']
                    transcript_id = annotations['transcript_id']
                    canon_transcript_ids = transcripts_by_gid.get(gene_id)
                    if not canon_transcript_ids:
                        genes_without_canon.add(gene_name)
                        continue
                    if not transcript_id in canon_transcript_ids:
                        continue

                start = int(start) - 1
                end = int(end)
                if end - start >= 3:
                    out.write('\t'.join([chrom, str(start), str(end), gene_name]) + '\n')
                    genes_set.add(gene_name)
                    region_cnt += 1

        if region_cnt % 10000 == 0:
            warn(f'Processed {len(genes_set)} genes, written {region_cnt} regions...')

    warn(f'Done. Processed {len(genes_set)} genes, written {region_cnt} regions')
    if genes_without_canon:
        warn(f'No canonical transcript for {len(genes_without_canon)} gene ids')


if __name__ == '__main__':
    main()

