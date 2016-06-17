#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

from collections import defaultdict, OrderedDict
import sys
from traceback import format_exc
from source.file_utils import adjust_path, verify_file, open_gzipsafe, add_suffix
from source.logger import err, info
from source.targetcov.Region import SortableByChrom
from source.utils import get_chr_lengths_from_seq, is_local

us_syn_path = '/ngs/reference_data/genomes/Hsapiens/common/HGNC_gene_synonyms.txt'
hg38_seq_fpath = '~/Dropbox/az/reference_data/hg38.fa' if is_local() else '/ngs/reference_data/genomes/Hsapiens/hg38/seq/hg38.fa'
hg19_seq_fpath = '~/Dropbox/az/reference_data/hg19.fa' if is_local() else '/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa'
canonical_hg38_transcripts_fpath = '~/Dropbox/az/reference_data/canonical_transcripts/canonical_transcripts_hg38.txt' if is_local() else '/ngs/reference_data/genomes/Hsapiens/hg38/canonical_transcripts.txt'
canonical_hg19_transcripts_fpath = '~/Dropbox/az/reference_data/canonical_transcripts/canonical_transcripts_hg19.txt' if is_local() else '/ngs/reference_data/genomes/Hsapiens/hg19/canonical_transcripts.txt'

ALL_EXONS = True
# MIR_AND_CDS_ONLY = False


'''
# * How to update reference data on servers *
# LOCAL:
make_exons.py hg19 RefSeq_knownGene.hg19.txt RefSeq.all_features.hg19.bed
make_exons.py hg38 RefSeq_knownGene.hg38.txt RefSeq.all_features.hg38.bed

egrep -w "miRNA|protein_coding" RefSeq.all_features.hg19.bed > RefSeq_CDS_miRNA.all_features.hg19.bed
egrep -w "miRNA|protein_coding" RefSeq.all_features.hg38.bed > RefSeq_CDS_miRNA.all_features.hg38.bed
egrep -w "miRNA|protein_coding" RefSeq.all_features.hg19.canon.bed > RefSeq_CDS_miRNA.all_features.hg19.canon.bed
egrep -w "miRNA|protein_coding" RefSeq.all_features.hg38.canon.bed > RefSeq_CDS_miRNA.all_features.hg38.canon.bed

grep -v "_hap" RefSeq.all_features.hg19.bed                 > RefSeq.all_features.hg19-noalt.bed
grep -v "_alt" RefSeq.all_features.hg38.bed                 > RefSeq.all_features.hg38-noalt.bed
grep -v "_hap" RefSeq.all_features.hg19.canon.bed           > RefSeq.all_features.hg19-noalt.canon.bed
grep -v "_alt" RefSeq.all_features.hg38.canon.bed           > RefSeq.all_features.hg38-noalt.canon.bed
grep -v "_hap" RefSeq_CDS_miRNA.all_features.hg19.bed       > RefSeq_CDS_miRNA.all_features.hg19-noalt.bed
grep -v "_alt" RefSeq_CDS_miRNA.all_features.hg38.bed       > RefSeq_CDS_miRNA.all_features.hg38-noalt.bed
grep -v "_hap" RefSeq_CDS_miRNA.all_features.hg19.canon.bed > RefSeq_CDS_miRNA.all_features.hg19-noalt.canon.bed
grep -v "_alt" RefSeq_CDS_miRNA.all_features.hg38.canon.bed > RefSeq_CDS_miRNA.all_features.hg38-noalt.canon.bed

grep -w CDS RefSeq_CDS_miRNA.all_features.hg19.canon.bed       | cut -f1,2,3,4 > RefSeq_CDS.hg19.bed
grep -w CDS RefSeq_CDS_miRNA.all_features.hg19-noalt.canon.bed | cut -f1,2,3,4 > RefSeq_CDS.hg19-noalt.bed
grep -w CDS RefSeq_CDS_miRNA.all_features.hg38.canon.bed       | cut -f1,2,3,4 > RefSeq_CDS.hg38.bed
grep -w CDS RefSeq_CDS_miRNA.all_features.hg38-noalt.canon.bed | cut -f1,2,3,4 > RefSeq_CDS.hg38-noalt.bed

grep -w chr20 RefSeq_CDS.hg19.bed > RefSeq_CDS.hg19-chr20.bed
grep -w chr20 RefSeq.all_features.hg19.bed > RefSeq.all_features.hg19-chr20.bed
grep -w chr20 RefSeq_CDS_miRNA.all_features.hg19.canon.bed > RefSeq_CDS_miRNA.all_features.hg19-chr20.canon.bed

# UK:
cp ~/Dropbox/az/reference_data/Exons/RefSeq/*.hg19*.bed /ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/RefSeq
cp ~/Dropbox/az/reference_data/Exons/RefSeq/*.hg38*.bed /ngs/reference_data/genomes/Hsapiens/hg38/bed/Exons/RefSeq
scp -r /ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/RefSeq/*.bed cniclhpc003:/ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/RefSeq
scp -r /ngs/reference_data/genomes/Hsapiens/hg38/bed/Exons/RefSeq/*.bed cniclhpc003:/ngs/reference_data/genomes/Hsapiens/hg38/bed/Exons/RefSeq

# US:
scp -r $uk:/ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/RefSeq/*.bed /ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/RefSeq
scp -r $uk:/ngs/reference_data/genomes/Hsapiens/hg38/bed/Exons/RefSeq/*.bed /ngs/reference_data/genomes/Hsapiens/hg38/bed/Exons/RefSeq
'''


def main():
    if len(sys.argv) < 4:
        info('The script writes all CDS, stop codon, and ncRNA exon regions for all known Ensembl genes, with associated gene symbols.')
        # info('When the gene name is found in HGNC, it get replaced with an approved name.         ')
        # info('If the gene is not charactirized (like LOC729737), this symbol is just kept as is.  ')
        info('                                                                                      ')
        info('Usage:                                                                                ')
        info('    ' + __file__ + ' hg19 db.gtf output.bed [HGNC_gene_synonyms.txt=' + us_syn_path + '] [additional_feature_list]')
        info('                                                                                      ')
        info('     where HGNC_gene_synonyms.txt (from http://www.genenames.org/cgi-bin/download) is:')
        info('     #Approved Symbol  Previous Symbols                    Synonyms                          Chromosome   Ensembl Gene ID   UCSC ID(supplied by UCSC)')
        info('     OR7E26P           OR7E67P, OR7E69P, OR7E70P, OR7E68P  OR1-51, OR1-72, OR1-73, OR912-95  19q13.43	    ENSG00000121410   uc002qsg.3')
        info('     ...                                                                              ')
        info('                                                                                      ')
        info('     or DB is Ensembl GTF ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz')
        info('     1  pseudogene            gene        11869  14412  .  +  .  gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene";')
        info('     1  processed_transcript  transcript  11869  14409  .  +  .  gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana";')
        info('     ...                                                                              ')
        info('                                                                                      ')
        info('     or DB is RefSeq GTF ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/H_sapiens/GFF/ref_GRCh38.p2_top_level.gff3.gz')
        info('     NC_000001.10    RefSeq          region       1       249250621       .       +       .       ID=id0;Name=1;Dbxref=taxon:9606;chromosome=1;gbkey=Src;genome=chromosome;mol_type=genomic DNA')
        info('     NC_000001.10    BestRefSeq      gene         11874   14409           .       +       .       ID=gene0;Name=DDX11L1;Dbxref=GeneID:100287102,HGNC:37102;description=DEAD%2FH %28Asp-Glu-Ala-Asp%2FHis%29 box helicase 11 like 1;gbkey=Gene;gene=DDX11L1;part=1%2F1;pseudo=true')
        info('     NC_000001.10    BestRefSeq      transcript   11874   14409           .       +       .       ID=rna0;Name=NR_046018.2;Parent=gene0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:37102;gbkey=misc_RNA;gene=DDX11L1;product=DEAD%2FH %28Asp-Glu-Ala-Asp%2FHis%29 box helicase 11 like 1;transcript_id=NR_046018.2')
        info('     NC_000001.10    BestRefSeq      exon         11874   12227           .       +       .       ID=id1;Parent=rna0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:37102;gbkey=misc_RNA;gene=DDX11L1;product=DEAD%2FH %28Asp-Glu-Ala-Asp%2FHis%29 box helicase 11 like 1;transcript_id=NR_046018.2')
        info('     ...                                                                              ')
        info('                                                                                      ')
        info('     or either RefSeq_knownGene.txt or UCSC_knownGene.txt (from http://genome.ucsc.edu/cgi-bin/hgTables) is:')
        info('     #hg19.knownGene.name  hg19.knownGene.chrom  hg19.knownGene.strand  hg19.knownGene.txStart  hg19.knownGene.txEnd  hg19.knownGene.exonCount  hg19.knownGene.exonStarts  hg19.knownGene.exonEnds  hg19.kgXref.geneSymbol')
        info('     uc001aaa.3	         chr1	               +	                  11873                   14409                 3                         11873,12612,13220,	      12227,12721,14409,	   DDX11L1')
        info('     ...                                                                              ')
        info('                                                                                      ')
        info('     Writes to Exons.bed                                                              ')
        info('                                                                                      ')
        info('See more info in http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Making+the+full+list+of+UCSC+exons+with+approved+HUGO+gene+symbols')
        sys.exit(1)

    genome_name = sys.argv[1]
    seq_fpath = hg19_seq_fpath if genome_name == 'hg19' else hg38_seq_fpath
    canonical_transcripts_fpath = canonical_hg19_transcripts_fpath if genome_name == 'hg19' else canonical_hg38_transcripts_fpath
    chr_lengths = get_chr_lengths_from_seq(seq_fpath)
    chr_order = {c: i for i, (c, l) in enumerate(chr_lengths)}

    input_fpath = verify_file(sys.argv[2])
    output_fpath = adjust_path(sys.argv[3])

    synonyms_fpath = None
    if len(sys.argv) > 4:
        synonyms_fpath = verify_file(sys.argv[4])
        info('Synonyms file provided ' + synonyms_fpath + '')
    else:
        info('No synonyms file provided, skipping approving')

    not_approved_fpath = None
    if len(sys.argv) > 5:
        not_approved_fpath = adjust_path(sys.argv[5])

    with open(verify_file(canonical_transcripts_fpath)) as f:
        canonical_transcripts_ids = set(l.strip().split('.')[0] for l in f)

    info('Reading the features...')
    with open_gzipsafe(input_fpath) as inp:
        l = inp.readline()
        if output_fpath.endswith('.gtf') or output_fpath.endswith('.gtf.gz'):
            gene_by_name_and_chrom = _proc_ensembl_gtf(inp, output_fpath, chr_order)
        elif output_fpath.endswith('.gff3') or output_fpath.endswith('.gff3.gz'):
            gene_by_name_and_chrom = _proc_refseq_gff3(inp, output_fpath, chr_order)
        else:
            gene_by_name_and_chrom = _proc_ucsc(inp, output_fpath, chr_order)

    if synonyms_fpath and synonyms_fpath != "''":
        gene_by_name_and_chrom, not_approved_gene_names = _approve(gene_by_name_and_chrom, synonyms_fpath)

        info('')
        info('Not approved by HGNC - ' + str(len(not_approved_gene_names)) + ' genes.')
        if not_approved_fpath:
            with open(not_approved_fpath, 'w') as f:
                f.write('#Searched as\tStatus\n')
                f.writelines((l + '\n' for l in not_approved_gene_names))
            info('Saved not approved to ' + not_approved_fpath)

        # with open('serialized_genes.txt', 'w') as f:
        #     for g in gene_by_name.values():
        #         f.write(str(g) + '\t' + str(g.db_id) + '\n')
        #         for e in g.exons:
        #             f.write('\t' + str(e) + '\n')

    info('Found:')
    info('  ' + str(len(gene_by_name_and_chrom)) + ' genes')

    genes = gene_by_name_and_chrom.values()

    coding_and_mirna_genes = [g for g in genes if any(t.biotype in ['protein_coding', 'miRNA'] for t in g.transcripts)]

    coding_genes = [g for g in coding_and_mirna_genes if any(t.biotype == 'protein_coding' for t in g.transcripts)]
    coding_transcripts = [t for g in coding_and_mirna_genes for t in g.transcripts if t.biotype == 'protein_coding']
    mirna_genes = [g for g in coding_and_mirna_genes if any(t.biotype == 'miRNA' for t in g.transcripts)]
    mirna_transcripts = [t for g in coding_and_mirna_genes for t in g.transcripts if t.biotype == 'miRNA']
    codingmiRNA_genes = [g for g in coding_and_mirna_genes if any(t.biotype == 'miRNA' for t in g.transcripts) and any(t.biotype == 'protein_coding' for t in g.transcripts)]
    info('  ' + str(len(coding_genes)) + ' coding genes')
    info('  ' + str(len(coding_transcripts)) + ' coding transcripts')
    info('  ' + str(len(mirna_genes)) + ' miRNA genes')
    info('  ' + str(len(mirna_transcripts)) + ' miRNA transcripts')
    info('  ' + str(len(codingmiRNA_genes)) + ' genes with both coding and miRNA transcripts')

    info()
    # info('Choosing genes with exons...')
    # genes = [g for g in genes if any(tx.exons for tx in g.transcripts)]
    # genes = [g for g in genes if any(tx.exons for tx in g.transcripts)]

    info('Choosing canonical...')
    canon_genes = choose_canonical(genes, canonical_transcripts_ids)

    info()
    info('Sorting and printing all regions...')
    print_genes(genes, output_fpath, canon_only=False)

    info()
    info('Sorting and printing canonical regions...')
    canon_output_fpath = add_suffix(output_fpath, 'canon')
    print_genes(canon_genes, canon_output_fpath, canon_only=True)

    info()
    info('Saved all regions to\n   ' + output_fpath + '\n   ' + canon_output_fpath)


def print_genes(genes, output_fpath, canon_only):
    regions = []
    already_added_gene_features = set()
    transcripts = []
    for g in genes:
        for tr in g.transcripts:
            if not canon_only or tr.is_canonical:
                transcripts.append(tr)
    for tr in sorted(transcripts, key=lambda _tr: _tr.get_key()):
        to_add_gene = all(tr2.biotype == 'protein_coding' for tr2 in tr.gene.transcripts if (tr2.is_canonical or not canon_only)) \
                      and tr.gene not in already_added_gene_features \
                      and (len(tr.gene.canonical_transcripts) == 1 or len(tr.gene.transcripts) == 1)
        if to_add_gene:
            # skip gene feature for all miRNA because there are multi-domain miRNA located in different
            # places with the same gene name
            regions.append(tr.gene)
            already_added_gene_features.add(tr.gene)
        if tr.exons:
            regions.append(tr)
            for e in tr.exons:
                regions.append(e)

    info('Writing ' + str(len(regions)) + ' regions')
    with open(adjust_path(output_fpath), 'w') as all_out:
        for r in regions:
            all_out.write(r.__str__())


def choose_canonical(genes, canonical_transcripts_ids):
    not_found_in_canon_coding_num = 0
    not_found_in_canon_coding_num_one_transcript = 0
    not_found_in_canon_mirna_num = 0
    not_found_in_canon_other_num = 0
    many_canon_coding_num = 0
    many_canon_mirna_num = 0
    many_canon_other_num = 0

    canon_genes = []
    for g in genes:
        _canon_tx = []
        for t in g.transcripts:
            if t.transcript_id in canonical_transcripts_ids:
                t.is_canonical = True
                _canon_tx.append(t)

        if len(_canon_tx) > 1:
            if any(t.biotype == 'protein_coding' for t in g.transcripts):
                many_canon_coding_num += 1
                longest_canon_t = max(_canon_tx, key=Transcript.length)
                for t in g.transcripts:
                    if t != longest_canon_t:
                        t.is_canonical = False
            elif any(t.biotype == 'miRNA' for t in g.transcripts):
                many_canon_mirna_num += 1
            else:
                many_canon_other_num += 1

        if len(_canon_tx) == 0:
            if any(t.biotype == 'protein_coding' for t in g.transcripts):
                not_found_in_canon_coding_num += 1
                if len(g.transcripts) == 1:
                    not_found_in_canon_coding_num_one_transcript += 1
                longest_t = max(g.transcripts, key=Transcript.length)
                longest_t.is_canonical = True
            elif any(t.biotype == 'miRNA' for t in g.transcripts):
                not_found_in_canon_mirna_num += 1
            else:
                not_found_in_canon_other_num += 1

        g.canonical_transcripts = [t for t in g.transcripts if t.is_canonical]
        if len(g.canonical_transcripts) > 0:
            if g.canonical_transcripts:
                canon_genes.append(g)

    info('Coding genes with canonical transcripts: ' +
         str(sum(1 for g in canon_genes if any(t.biotype == 'protein_coding' for t in g.canonical_transcripts))))
    info('Coding canonical transcripts: ' +
         str(sum(1 for g in canon_genes for t in g.canonical_transcripts if t.biotype == 'protein_coding')))
    info('miRNA genes with canonical transcripts: ' +
         str(sum(1 for g in canon_genes if any(t.biotype == 'miRNA' for t in g.canonical_transcripts))))
    info('miRNA canonical transcripts: ' +
         str(sum(1 for g in canon_genes for t in g.canonical_transcripts if t.biotype == 'miRNA')))
    info('Other canonical transcripts: ' +
         str(sum(1 for g in canon_genes for t in g.canonical_transcripts if t.biotype not in ['miRNA', 'protein_coding'])))
    info('Other genes with canonical transcripts: ' +
         str(sum(1 for g in canon_genes if any(t.biotype not in ['miRNA', 'protein_coding'] for t in g.canonical_transcripts))))

    info()
    info('Coding genes with no canonical transcripts (picking longest out of the rest): ' + str(not_found_in_canon_coding_num))
    info('miRNA genes with no canonical transcripts (skipping all): ' + str(not_found_in_canon_mirna_num))
    info('Other genes with no canonical transcripts (skipping all): ' + str(not_found_in_canon_other_num))
    info('Coding genes with many canonical transcripts (picking longest): ' + str(many_canon_coding_num))
    info('miRNA genes with many canonical transcripts (keeping all): ' + str(many_canon_mirna_num))
    info('Other genes with many canonical transcripts (keeping all): ' + str(many_canon_other_num))

    return canon_genes


def _approve(gene_by_name, synonyms_fpath):
    approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym = \
        read_approved_genes(synonyms_fpath)

    not_approved_gene_names = list()
    gene_after_approving_by_name = OrderedDict()
    total_approved = 0
    total_not_approved = 0
    j = 0
    for g in gene_by_name.values():
        if len(g.exons) == 0:
            continue

        gene_after_approving_by_name[g.name] = g
        if is_approved_symbol(g.name, approved_gene_by_name):
            gene_after_approving_by_name[g.name] = g
            total_approved += 1
        else:
            not_approved_gene_names.append(g.name)
            total_not_approved += 1

        j += 1
        if j % 1000 == 0:
            info('processed ' + str(j / 1000) + 'k genes...')

    info('-----')
    info('Total: ' + str(j))
    if approved_gene_by_name:
        info('Total approved: ' + str(total_approved))
        info('Total not approved: ' + str(total_not_approved))
    info()
    info('Saving genes...')

    gene_features = 0
    features_counter = defaultdict(int)
    biotypes_counter = defaultdict(int)
    no_exon_gene_num = 0

    filtered_gene_after_approving_by_name = OrderedDict()
    for g in gene_after_approving_by_name.values():
        if len(g.exons) == 0:
            no_exon_gene_num += 1
        else:
            filtered_gene_after_approving_by_name[g.name] = g

            gene_features += 1
            features_counter[g.feature] += 1
            biotypes_counter[g.biotype] += 1

            for e in g.exons:
                features_counter[e.feature] += 1

                if e.feature == 'exon': e.feature = 'Exon'
                elif e.feature == 'stop_codon': e.feature = 'CDS'
                else: e.feature = e.feature[0].upper() + e.feature[1:]

    info('Skipped {} genes with no sub-features.'.format(no_exon_gene_num))
    info('Approved {} genes, including:'.format(gene_features))
    info('    Gene: {}'.format(features_counter['Gene']))
    info('    Multi_Gene: {}'.format(features_counter['Multi_Gene']))
    info('')

    info('Out of total: {} protein coding genes, {} ncRNA genes, including:'.format(
        biotypes_counter['protein_coding'], sum(biotypes_counter.values()) - biotypes_counter['protein_coding']))
    for bt, cnt in biotypes_counter.items():
        if bt != 'protein_coding':
            err('    ' + bt + ': ' + str(cnt))

    info()
    if ALL_EXONS:
        info('Found {} exons.'.format(features_counter['exon']))
    else:
        info('Also found {} CDS, {} stop codons, and {} ncRNA exons.'.format(
            features_counter['CDS'], features_counter['stop_codon'], features_counter['exon']))

    return filtered_gene_after_approving_by_name, not_approved_gene_names


class ApprovedGene:
    def __init__(self, name, prev_names, synonyms, chrom, ucsc_id=None, ensembl_id=None):
        self.name = name
        self.prev_names = prev_names
        self.synonyms = synonyms
        self.chrom = chrom

        self.db_id = ensembl_id


def parse_hgnc_chrom(chrom):
    if chrom in ['reserved', 'c10_B']:
        return None

    CHROMS = ['Y', 'X', 'mitochondria']
    for i in range(22, 0, -1):
        CHROMS.append(str(i))

    for c in CHROMS:
        if chrom.startswith(c):
            if c == 'mitochondria':
                return 'chrM'
            return 'chr' + c

    info('  Notice: cannot parse chromosome ' + chrom)
    return None


def parse_ensembl_chrom(chrom):
    CHROMS = ['Y', 'X', 'MT']
    for i in range(22, 0, -1):
        CHROMS.append(str(i))

    for c in CHROMS:
        if chrom.startswith(c):
            if c == 'MT':
                return 'chrM'
            return 'chr' + c

    return None


def read_approved_genes(synonyms_fpath):
    approved_gene_by_name = dict()
    approved_gnames_by_prev_gname = defaultdict(list)
    approved_gnames_by_synonym = defaultdict(list)

    info('Parsing HGNC database ' + synonyms_fpath + '...')
    with open(synonyms_fpath) as f:
        i = 0
        for l in f:
            if l and not l.startswith('#'):
                approved_gn, prev_names, synonyms, hgnc_chrom, ensembl_id, ucsc_id = l.replace('\n', '').split('\t')
                if hgnc_chrom:
                    hgnc_chrom = parse_hgnc_chrom(hgnc_chrom)

                approved_gene = ApprovedGene(approved_gn, prev_names, synonyms, hgnc_chrom, ucsc_id, ensembl_id)
                approved_gene_by_name[approved_gn] = approved_gene

                for gn in prev_names.split(', '):
                    if gn:
                        approved_gnames_by_prev_gname[gn].append(approved_gene)

                for gn in synonyms.split(', '):
                    if gn:
                        approved_gnames_by_synonym[gn].append(approved_gene)
            i += 1
        info('  Processed ' + str(i) + ' lines from ' + synonyms_fpath)
        info()

    return approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym


def _check_gene_symbol(approved_gene, gene_symbol, db_id, chrom):
    if db_id and db_id != approved_gene.db_id:
        # sys.stderr.write('Discordant db ids for ' + gene_symbol + ': db id = ' +
        # str(db_id) + ', in HGNC it is ' + str(approved_gene.db_id) + '\n')
        pass
    # else:
        # sys.stderr.write('Accordant db ids for ' + gene_symbol + ': db id = ' +
        # str(db_id) + ', in HGNC it is ' + str(approved_gene.db_id) + '\n')

    if chrom and chrom != approved_gene.chrom:
        # sys.stderr.write('Discordant chroms for ' + gene_symbol + ': chrom = ' +
        # chrom + ', in HGNC chrom is ' + approved_gene.chrom + '\n')
        return None

    return approved_gene


def get_approved_gene_symbol(approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym,
                             gene_symbol, db_id='', db_chrom='', indent=''):
    if gene_symbol in approved_gene_by_name:
        if _check_gene_symbol(approved_gene_by_name[gene_symbol], gene_symbol, db_id, db_chrom):
            return approved_gene_by_name[gene_symbol].name, None

    info(indent + 'Gene name ' + gene_symbol + ' is not approved, searching for an approved version... ',
        ending='', print_date=False)

    def _get_approved_genes_by_kind(approved_genes, kind):
        if not approved_genes:
            return 'NOT FOUND'

        if len(approved_genes) > 1:
            approved_genes_same_ucsc = [g for g in approved_genes if g.db_id == db_id]

            if len(approved_genes_same_ucsc) > 1:
                err(' ERROR: multiple approved gene names for ' + gene_symbol + ' (as ' + kind + ') with ucsc_id ' +
                    db_id + ': ' + ', '.join(g.name for g in approved_genes_same_ucsc) + '', print_date=False)
                return 'AMBIGUOUS'

            if len(approved_genes_same_ucsc) == 1:
                if _check_gene_symbol(approved_genes_same_ucsc[0], gene_symbol, db_id, db_chrom):
                    err(' found approved gene for ' + gene_symbol + ' (as ' + kind + ') with ucsc_id ' + db_id,
                        print_date=False)
                    return approved_genes_same_ucsc[0].name

            # Ok, no genes with same ucsc id, or not the same chromosome for them.

            approved_genes_same_chrom = [g for g in approved_genes if g.chrom == db_chrom]

            if len(approved_genes_same_chrom) > 1:
                err(' ERROR: multiple approved gene names for ' + gene_symbol + ' (as ' + kind + ') with chrom ' +
                    db_chrom + ', '.join(g.name for g in approved_genes_same_ucsc) + '', print_date=False)
                return 'AMBIGUOUS'

            if len(approved_genes_same_chrom) == 1:
                g = approved_genes_same_chrom[0]
                info(' only ' + g.name + ' for ' + gene_symbol + ' (as ' + kind + ') has the same chrom '
                    + db_chrom + ', picking it', print_date=False)
                if _check_gene_symbol(g, gene_symbol, db_id, db_chrom):
                    return g.name
                else:
                    return 'NOT FOUND'

            if len(approved_genes_same_chrom) == 0:
                err(' ERROR: no approved gene names for ' + gene_symbol + ' (as ' + kind + ') with same chrom '
                    + db_chrom + '', print_date=False)
                return 'NOT FOUND'

        if len(approved_genes) == 1:
            if _check_gene_symbol(approved_genes[0], gene_symbol, db_id, db_chrom):
                info(' found approved gene symbol for ' + gene_symbol + ': ' + approved_genes[0].name + ' (as '
                    + kind + ')', print_date=False)
                return approved_genes[0].name

        return 'NOT FOUND'

    res = _get_approved_genes_by_kind(approved_gnames_by_prev_gname.get(gene_symbol), 'prev')
    if res == 'AMBIGUOUS':
        return None, 'AMBIGUOUS\tAS PREV'
    elif res == 'NOT FOUND':
        res = _get_approved_genes_by_kind(approved_gnames_by_synonym.get(gene_symbol), 'synonym')
        if res == 'AMBIGUOUS':
            return None, res + '\tAS SYNONYM'
        if res == 'NOT FOUND':
            err(' not found.', print_date=False)
            return None, res
        else:
            info(indent + 'Finally found approved gene for ' + gene_symbol + ' (as synonym): ' + res, print_date=False)
            return res, None
    else:
        info(indent + 'Finally found approved gene for ' + gene_symbol + ' (as prev): ' + res, print_date=False)
        return res, None


def _proc_ucsc(inp, output_fpath, chr_order):  #, approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym):
    gene_by_name_and_chrom = dict()

    prev_chrom = None
    for l in inp:
        if l and not l.startswith('#'):
            transcript_id, ucsc_chrom, strand, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, gene_symbol =\
                l.replace('\n', '').split('\t')
            cdsStart = int(cdsStart)
            cdsEnd = int(cdsEnd)
            exonCount = int(exonCount)
            exonStarts = [int(v) + 1 for v in exonStarts.split(',') if v]
            exonEnds = map(int, filter(None, exonEnds.split(',')))

            if ucsc_chrom != prev_chrom:
                info(ucsc_chrom)
                prev_chrom = ucsc_chrom

            # approved_gene_symbol, status = get_approved_gene_symbol(
            #     approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym,
            #     gene_symbol, ucsc_id, ucsc_chrom)
            #
            # if not approved_gene_symbol:
            #     not_approved_gene_names.append(gene_symbol + '\t' + status)
            #     if DO_APPROVE:
            #         continue
            #     else:
            #         approved_gene_symbol = gene_symbol

            txStart = exonStarts[0] - 1
            txEnd = exonEnds[exonCount - 1]

            # out.write('\t'.join([ucsc_chrom, str(min(txStart, cdsStart)), str(max(txEnd, cdsEnd)),
            #                      gene_symbol, '.', strand, 'Gene', '.']) + '\n')

            assert txStart <= cdsStart, l
            assert txEnd >= cdsEnd, l

            if (gene_symbol, ucsc_chrom) not in gene_by_name_and_chrom:
                gene = Gene(ucsc_chrom, chr_order.get(ucsc_chrom), gene_symbol, strand)
                gene_by_name_and_chrom[(gene_symbol, ucsc_chrom)] = gene
            gene = gene_by_name_and_chrom[(gene_symbol, ucsc_chrom)]

            transcript = Transcript(gene, transcript_id, txStart, txEnd, strand)
            gene.transcripts.append(transcript)

            if gene_symbol == 'AKT1':
                pass

            for exon_number, eStart, eEnd in zip(
                   range(exonCount),
                   [s for s in exonStarts if s],
                   [e for e in exonEnds if e]):
                eStart -= 1

                exon = None
                if eEnd <= cdsStart or eStart > cdsEnd:  # usually it means cdsStart = 0,
                                                         # no CDS for this gene, thus reporting exons
                    if gene_symbol.startswith('MIR'):
                        if transcript.biotype != 'protein_coding':
                            transcript.biotype = 'miRNA'
                    else:
                        if transcript.biotype not in ['protein_coding', 'miRNA']:
                            transcript.biotype = 'other'
                    exon = Exon(transcript, eStart, eEnd, 'Exon', exon_number)

                else:
                    transcript.biotype = 'protein_coding'
                    if cdsStart <= eStart:
                        exon = Exon(transcript, eStart, eEnd, 'CDS', exon_number)
                    elif eEnd > cdsStart:
                        exon = Exon(transcript, cdsStart, eEnd, 'CDS', exon_number)
                    else:
                        err('Warn: exon ' + str(eStart) + ':' + str(eEnd) + ' does not contain CDS, CDS start = ' + str(cdsStart))

                if exon:
                    transcript.exons.append(exon)

    return gene_by_name_and_chrom


class Transcript(SortableByChrom):
    def __init__(self, gene, transcript_id, start, end, strand, biotype=None):
        SortableByChrom.__init__(self, gene.chrom, gene.chrom_ref_order)
        self.gene = gene
        self.transcript_id = transcript_id
        self.start = start
        self.end = end
        self.strand = strand
        self.biotype = biotype
        self.feature = 'Transcript'

        self.exons = []
        self.is_canonical = False

    def __str__(self):
        fs = [self.chrom,
              '{}'.format(self.start) if self.start is not None else '.',
              '{}'.format(self.end) if self.end is not None else '.',
              self.gene.name or '.',
              '.',
              self.strand or '.',
              self.feature or '.',
              self.biotype or '.',
              self.transcript_id,
              ]
        return '\t'.join(fs) + '\n'

    def get_key(self):
        return self.chrom_ref_order, self.start, self.end

    def length(self):
        return sum(e.end - e.start for e in self.exons)


class Gene(SortableByChrom):
    def __init__(self, chrom, chrom_ref_order, name, strand, biotype='', db_id='', source=''):
        SortableByChrom.__init__(self, chrom, chrom_ref_order)
        self.name = name
        self.strand = strand
        self.biotype = biotype
        self.db_id = db_id
        self.feature = 'Gene'
        self.source = source

        self.approved_gname = None

        self.transcripts = []
        self.canonical_transcripts = []

    def get_start(self):
        if len(self.canonical_transcripts) == 1:
            return self.canonical_transcripts[0].start
        elif len(self.transcripts) == 1:
            return self.transcripts[0].start
        else:
            return 'There must be exactly 1 transcript or exactly 1 canonical transcript in a gene to get gene start. ' \
                   'Number of transcripts is ' + str(len(self.transcripts)) + ', ' \
                   'Number of canonical transcripts is ' + str(len(self.canonical_transcripts))

    def get_end(self):
        if len(self.canonical_transcripts) == 1:
            return self.canonical_transcripts[0].end
        elif len(self.transcripts) == 1:
            return self.transcripts[0].end
        else:
            return 'There must be exactly 1 transcript or exactly 1 canonical transcript in a gene to get gene start. ' \
                   'Number of transcripts is ' + str(len(self.transcripts)) + ', ' \
                   'Number of canonical transcripts is ' + str(len(self.canonical_transcripts))

    def __str__(self):
        fs = [self.chrom,
              '{}'.format(self.get_start()) if self.get_start() is not None else '.',
              '{}'.format(self.get_end()) if self.get_end() is not None else '.',
              self.name or '.',
              '.',
              self.strand or '.',
              self.feature or '.',
              self.biotype or '.',
              '.',
              ]
        return '\t'.join(fs) + '\n'

    def __repr__(self):
        return '{self.name} {self.chrom}:{self.start}-{self.end} {self.biotype} ' \
               '{self.db_id} {self.source}'.format(self=self)

    def get_key(self):
        return self.chrom_ref_order, self.get_start(), self.get_end()


class Exon(SortableByChrom):
    def __init__(self, transcript, start, end, feature, exon_number):
        SortableByChrom.__init__(self, transcript.gene.chrom, transcript.gene.chrom_ref_order)
        self.transcript = transcript
        self.start = start
        self.end = end
        self.feature = feature
        self.exon_number = exon_number

    def __str__(self):
        fs = [self.transcript.gene.chrom,
              '{}'.format(self.start) if self.start is not None else '.',
              '{}'.format(self.end) if self.end is not None else '.',
              self.transcript.gene.name or '.',
              '{}'.format(self.exon_number) if self.exon_number is not None else '.',
              self.transcript.gene.strand or '.',
              self.feature or '.',
              self.transcript.biotype or '.',
              self.transcript.transcript_id or '.',
             ]
        return '\t'.join(fs) + '\n'

    def get_key(self):
        return self.chrom_ref_order, self.start, self.end


def is_approved_symbol(gname, approved_gene_by_name):
    if gname not in approved_gene_by_name:
        # gname2 = gname.split('.')[0]
        # if gname != gname2:
        #     if gname2 not in approved_gene_by_name:
        return False
    return True


def _proc_refseq_gff3(inp, out, chr_order, additional_feature_list=None):
    gene_by_name = OrderedDict()
    gene_by_id = OrderedDict()

    info('Parsing RefSeq GFF3...')
    total_lines = 0
    total_non_coding_genes = 0

    # TODO: from here

    for l in inp:
        if l and not l.startswith('#'):
            # refseq_id, db, feature, start, end,

            chrom, _, feature, start, end, _, strand, _, props_line = l.replace('\n', '').split('\t')

            # if is_local():
            #     if chrom != '21':
            #         continue

            total_lines += 1
            if total_lines % 1000 == 0:
                info(str(total_lines / 1000) + 'k lines, ' + str(len(gene_by_name)) + ' genes found')
                sys.stdout.flush()

            try:
                _prop_dict = dict((t.strip().split(' ')[0], ' '.join(t.strip().split(' ')[1:]))
                                  for t in props_line.split(';') if t.strip())
            except ValueError:
                sys.stderr.write(format_exc())
                sys.stderr.write(l)

            gene_symbol = _rm_quotes(_prop_dict['gene_name'])
            gene_id = _rm_quotes(_prop_dict['gene_id'])
            gene_biotype = _rm_quotes(_prop_dict['gene_biotype'])
            gene_source = _rm_quotes(_prop_dict['gene_source'])

            # if gene_symbol == 'PTENP1':
            #     sys.stderr.write('PTENP1\n')

            if not ALL_EXONS and gene_biotype not in [
                'protein_coding',
                'nonsense_mediated_decay',
                'non_stop_decay',
                'processed_transcript',
                'polymorphic_pseudogene',
                'sense_intronic',
                'sense_overlapping',
                'antisense',

            ] and not any(b in gene_biotype for b in ['RNA', 'IG_', 'TR_']):
                total_non_coding_genes += 1
                continue

            full_feature_list = ['gene', 'CDS', 'stop_codon', 'exon'] + additional_feature_list
            if ALL_EXONS:
                full_feature_list = ['gene', 'exon']
            # sys.stderr.write('Full feature list: ' + str(full_feature_list) + '\n')
            if feature not in full_feature_list:
                continue

            start, end = int(start) - 1, int(end)

            if int(end) <= int(start):
                info('Error: start > end: ' + l)
                continue

            chrom = parse_ensembl_chrom(chrom)
            if not chrom:
                continue

            if feature == 'gene':
                # assert gene_biotype == biotype, 'Gene: gene_biotype "' + gene_biotype + '"
                # do not match biotype "' + biotype + '" for ' + gene_symbol

                gene = Gene(chrom, chr_order.get(chrom), start, end, gene_symbol, strand,
                            gene_biotype, gene_id, gene_source)

                if gene.name in gene_by_name:
                    prev_gene = gene_by_name[gene.name]

                    if gene.source != prev_gene.source:
                        err('    Duplicated gene in different databases:')
                        err('        This: ' + gene.__repr__())
                        err('        Prev: ' + prev_gene.__repr__())
                        # answer = raw_input('Which one to pick? This (1), prev (2), longest (Enter): ')
                        #
                        # if answer == '1' or answer == '' and gene.end - gene.start >
                        # prev_gene.end - prev_gene.start:
                        #     del gene_by_name[prev_gene.name]
                        #     del gene_by_id[prev_gene.db_id]
                        #
                        # else:
                        #     continue

                        if gene.source == 'ensembl' or prev_gene.source == 'havana':
                            del gene_by_name[prev_gene.name]
                            del gene_by_id[prev_gene.db_id]
                            err('        Picking up this one.')

                        if prev_gene.source == 'ensembl' or gene.source == 'havana':
                            err('        Picking up previous one.')
                            continue

                    else:
                        err('    Duplicated gene in ' + gene.source + ':')
                        err('        ' + gene.__repr__())
                        prev_gene.start = min(prev_gene.start, gene.start)
                        prev_gene.end = max(prev_gene.end, gene.end)
                        prev_gene.feature = 'Multi_Gene'
                        continue

                    err('')

                gene_by_name[gene_symbol] = gene
                gene_by_id[gene_id] = gene

            elif feature in ['CDS', 'stop_codon'] \
                    or feature == 'exon' and ('RNA' in gene_biotype or ALL_EXONS) \
                    or feature in additional_feature_list:
                assert gene_symbol in gene_by_name, 'Error: ' + feature + ' record before gene record ' + \
                        gene_symbol + ', ' + gene_id + '; gene_by_name: ' + str(gene_by_name.keys())
                gene = gene_by_name[gene_symbol]
                if gene.gene_id == gene_id:
                    assert gene_biotype == gene.biotype, feature + ': gene_biotype "' + gene_biotype + \
                         '" do not match biotype "' + gene.biotype + '" for ' + gene_symbol
                    exon = Exon(gene, start, end, gene_biotype, feature)
                    gene.exons.append(exon)

    info()
    info(
        'Processed ' +
        str(total_lines) + ' lines, ' +
        str(total_non_coding_genes) + ' non-coding genes skipped, ' +
        str(len(gene_by_name)) + ' coding genes found')
    info()
    return gene_by_name


def _proc_ensembl_gtf(inp, out, chr_order, additional_feature_list=None):
    if additional_feature_list is None:
        additional_feature_list = []

    info('additional_feature_list = ' + str(additional_feature_list))

    gene_by_name = OrderedDict()
    gene_by_id = OrderedDict()

    info('Parsing Ensembl input...')
    total_lines = 0
    total_non_coding_genes = 0

    for l in inp:
        if l and not l.startswith('#'):
            chrom, _, feature, start, end, _, strand, _, props_line = l.replace('\n', '').split('\t')

            # if is_local():
            #     if chrom != '21':
            #         continue

            total_lines += 1
            if total_lines % 1000 == 0:
                info(str(total_lines / 1000) + 'k lines, ' + str(len(gene_by_name)) + ' genes found')
                sys.stdout.flush()

            try:
                _prop_dict = dict((t.strip().split(' ')[0], ' '.join(t.strip().split(' ')[1:]))
                                  for t in props_line.split(';') if t.strip())
            except ValueError:
                sys.stderr.write(format_exc())
                sys.stderr.write(l)

            gene_symbol = _rm_quotes(_prop_dict['gene_name'])
            gene_id = _rm_quotes(_prop_dict['gene_id'])
            gene_biotype = _rm_quotes(_prop_dict['gene_biotype'])
            gene_source = _rm_quotes(_prop_dict['gene_source'])

            # if gene_symbol == 'PTENP1':
            #     sys.stderr.write('PTENP1\n')

            if not ALL_EXONS and gene_biotype not in [
                'protein_coding',
                'nonsense_mediated_decay',
                'non_stop_decay',
                'processed_transcript',
                'polymorphic_pseudogene',
                'sense_intronic',
                'sense_overlapping',
                'antisense',

            ] and not any(b in gene_biotype for b in ['RNA', 'IG_', 'TR_']):
                total_non_coding_genes += 1
                continue

            full_feature_list = ['gene', 'CDS', 'stop_codon', 'exon'] + additional_feature_list
            if ALL_EXONS:
                full_feature_list = ['gene', 'exon']
            # sys.stderr.write('Full feature list: ' + str(full_feature_list) + '\n')
            if feature not in full_feature_list:
                continue

            start, end = int(start) - 1, int(end)

            if int(end) <= int(start):
                info('Error: start > end: ' + l)
                continue

            chrom = parse_ensembl_chrom(chrom)
            if not chrom:
                continue

            if feature == 'gene':
                # assert gene_biotype == biotype, 'Gene: gene_biotype "' + gene_biotype + '"
                # do not match biotype "' + biotype + '" for ' + gene_symbol

                gene = Gene(chrom, chr_order.get(chrom), start, end, gene_symbol, strand,
                            gene_biotype, gene_id, gene_source)

                if gene.name in gene_by_name:
                    prev_gene = gene_by_name[gene.name]

                    if gene.source != prev_gene.source:
                        err('    Duplicated gene in different databases:')
                        err('        This: ' + gene.__repr__())
                        err('        Prev: ' + prev_gene.__repr__())
                        # answer = raw_input('Which one to pick? This (1), prev (2), longest (Enter): ')
                        #
                        # if answer == '1' or answer == '' and gene.end - gene.start >
                        # prev_gene.end - prev_gene.start:
                        #     del gene_by_name[prev_gene.name]
                        #     del gene_by_id[prev_gene.db_id]
                        #
                        # else:
                        #     continue

                        if gene.source == 'ensembl' or prev_gene.source == 'havana':
                            del gene_by_name[prev_gene.name]
                            del gene_by_id[prev_gene.db_id]
                            err('        Picking up this one.')

                        if prev_gene.source == 'ensembl' or gene.source == 'havana':
                            err('        Picking up previous one.')
                            continue

                    else:
                        err('    Duplicated gene in ' + gene.source + ':')
                        err('        ' + gene.__repr__())
                        prev_gene.start = min(prev_gene.start, gene.start)
                        prev_gene.end = max(prev_gene.end, gene.end)
                        prev_gene.feature = 'Multi_Gene'
                        continue

                    err('')

                gene_by_name[gene_symbol] = gene
                gene_by_id[gene_id] = gene

            elif feature in ['CDS', 'stop_codon'] \
                    or feature == 'exon' and ('RNA' in gene_biotype or ALL_EXONS) \
                    or feature in additional_feature_list:
                assert gene_symbol in gene_by_name, 'Error: ' + feature + ' record before gene record ' + \
                        gene_symbol + ', ' + gene_id + '; gene_by_name: ' + str(gene_by_name.keys())
                gene = gene_by_name[gene_symbol]
                if gene.gene_id == gene_id:
                    assert gene_biotype == gene.biotype, feature + ': gene_biotype "' + gene_biotype + \
                         '" do not match biotype "' + gene.biotype + '" for ' + gene_symbol
                    exon = Exon(gene, start, end, gene_biotype, feature)
                    gene.exons.append(exon)

    info()
    info(
        'Processed ' +
        str(total_lines) + ' lines, ' +
        str(total_non_coding_genes) + ' non-coding genes skipped, ' +
        str(len(gene_by_name)) + ' coding genes found')
    info()
    return gene_by_name


if __name__ == '__main__':
    main()
