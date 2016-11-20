#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

from collections import defaultdict, OrderedDict
import sys
from traceback import format_exc
from source.file_utils import adjust_path

DO_APPROVE = False
ALL_EXONS = True


def main():
    if len(sys.argv) < 2:
        sys.stderr.write('The script writes all CDS, stop codon, and ncRNA exon regions for all known Ensembl genes, with associated gene symbols.\n')
        sys.stderr.write('When the gene name is found in HGNC, it get replaced with an approved name.\n')
        sys.stderr.write('If the gene is not charactirized (like LOC729737), this symbol is just kept as is.\n')
        sys.stderr.write('\n')
        sys.stderr.write('Usage:\n')
        sys.stderr.write('    ' + __file__ + ' Ensembl.gtf [HGNC_cBio_genes.tsv] [additional_feature_list] > Exons.bed\n')
        sys.stderr.write('\n')
        sys.stderr.write('   where HGNC_gene_synonyms.txt (from http://www.genenames.org/cgi-bin/download) is:\n')
        sys.stderr.write('     #Approved Symbol  Previous Symbols                    Synonyms                          Chromosome   Ensembl Gene ID   UCSC ID(supplied by UCSC)\n')
        sys.stderr.write('     OR7E26P           OR7E67P, OR7E69P, OR7E70P, OR7E68P  OR1-51, OR1-72, OR1-73, OR912-95  19q13.43	    ENSG00000121410   uc002qsg.3\n')
        sys.stderr.write('     ...\n')
        sys.stderr.write('\n')
        sys.stderr.write('   feature_list is by default empty, but could be transcript\n')
        sys.stderr.write('\n')
        sys.stderr.write('   and UCSC_knownGene.txt (from http://genome.ucsc.edu/cgi-bin/hgTables) is:\n')
        sys.stderr.write('     #hg19.knownGene.name  hg19.knownGene.chrom  hg19.knownGene.strand  hg19.knownGene.txStart  hg19.knownGene.txEnd  hg19.knownGene.exonCount  hg19.knownGene.exonStarts  hg19.knownGene.exonEnds  hg19.kgXref.geneSymbol\n')
        sys.stderr.write('     uc001aaa.3	          chr1	                +	                   11873                   14409                 3                         11873,12612,13220,	      12227,12721,14409,	   DDX11L1\n')
        sys.stderr.write('     ...\n')
        sys.stderr.write('   or Ensembl.gtf (ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz)')
        sys.stderr.write('     1  pseudogene            gene        11869  14412  .  +  .  gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene";')
        sys.stderr.write('     1  processed_transcript  transcript  11869  14409  .  +  .  gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana";')
        sys.stderr.write('     ...\n')
        sys.stderr.write('\n')
        sys.stderr.write('   Writes to Exons.bed\n')
        sys.stderr.write('\n')
        sys.stderr.write('See more info in http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Making+the+full+list+of+UCSC+exons+with+approved+HUGO+gene+symbols\n')
        sys.exit(1)

    # if is_local():
    #     sys.stderr.write('Local: will run only for chr21\n')
    #     sys.stderr.write('\n')

    input_fpath = adjust_path(sys.argv[1])
    hgnc_fpath = adjust_path(sys.argv[2])
    approved_gene_by_name = None
    if hgnc_fpath and hgnc_fpath != "''":
        sys.stderr.write('Synonyms file provided ' + hgnc_fpath + '\n')
        approved_gene_by_name = read_hgnc_genes(hgnc_fpath)
    else:
        sys.stderr.write('No synonyms file provided, skipping approving\n')

    out = sys.stdout
    with open(input_fpath) as inp:
        _ = inp.readline()
        not_approved_gene_names = _proc_ensembl(inp, out, approved_gene_by_name)

    # sys.stderr.write('\n')
    # sys.stderr.write('Not approved by HGNC - ' + str(len(not_approved_gene_names.keys())) + ' genes.\n')
    # if not_approved_fpath:
    #     with open(not_approved_fpath, 'w') as f:
    #         f.write('#Searched as\tStatus\n')
    #         f.writelines((gn + '\t' + st + '\n' for gn, st in not_approved_gene_names.items()))
    #     sys.stderr.write('Saved not approved to ' + not_approved_fpath + '\n')


class ApprovedGene:
    def __init__(self, symbol, cytoband, locus_group, entrez_gene_id, ensembl_id):
        self.name = symbol
        self.cytoband = cytoband
        self.entrez_gene_id = entrez_gene_id
        self.type = locus_group
        self.ensembl_id = ensembl_id
        self.length = None
        self.source = None

    def __str__(self):
        return '\t'.join([
            self.entrez_gene_id,
            self.name,
            self.type,
            self.cytoband,
            str(self.length)
        ])


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

    sys.stderr.write('  Cannot parse chromosome ' + chrom + '\n')
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


def read_hgnc_genes(hgnc_fpath):
    approved_gene_by_name = dict()

    sys.stderr.write('Parsing HGNC database ' + hgnc_fpath + '...\n')
    with open(hgnc_fpath) as f:
        i = 0
        for l in f:
            if l and not l.startswith('#'):
                symbol, locus_group, cytoband, entrez_gene_id, ensembl_id, entrez_gene_id_ncbi = l.replace('\n', '').split('\t')
                approved_gene = ApprovedGene(symbol, cytoband, locus_group, entrez_gene_id_ncbi, ensembl_id)
                approved_gene_by_name[symbol] = approved_gene

            i += 1
        sys.stderr.write('  Processed ' + str(i) + ' lines from ' + hgnc_fpath + '\n')
        sys.stderr.write('\n')

    return approved_gene_by_name


def _check_gene_symbol(approved_gene, gene_symbol, db_id, chrom):
    if db_id and db_id != approved_gene.db_id:
        # sys.stderr.write('Discordant db ids for ' + gene_symbol + ': db id = ' + str(db_id) + ', in HGNC it is ' + str(approved_gene.db_id) + '\n')
        pass
    # else:
        # sys.stderr.write('Accordant db ids for ' + gene_symbol + ': db id = ' + str(db_id) + ', in HGNC it is ' + str(approved_gene.db_id) + '\n')

    if chrom and chrom != approved_gene.chrom:
        # sys.stderr.write('Discordant chroms for ' + gene_symbol + ': chrom = ' + chrom + ', in HGNC chrom is ' + approved_gene.chrom + '\n')
        return None

    return approved_gene


def get_approved_gene_symbol(approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym,
                              gene_symbol, db_id='', db_chrom=''):
    if gene_symbol in approved_gene_by_name:
        if _check_gene_symbol(approved_gene_by_name[gene_symbol], gene_symbol, db_id, db_chrom):
            return approved_gene_by_name[gene_symbol].name, None

    sys.stderr.write('Gene name ' + gene_symbol + ' is not approved, searching for an approved version.\n')

    def _get_approved_genes_by_kind(approved_genes, kind):
        if not approved_genes:
            return 'NOT FOUND'

        if len(approved_genes) > 1:
            approved_genes_same_ucsc = [g for g in approved_genes if g.db_id == db_id]

            if len(approved_genes_same_ucsc) > 1:
                sys.stderr.write('  Error: multiple approved gene names for ' + gene_symbol + ' (as ' + kind + ') with ucsc_id ' + db_id + ': ' + ', '.join(g.name for g in approved_genes_same_ucsc) + '\n')
                return 'AMBIGUOUS'

            if len(approved_genes_same_ucsc) == 1:
                if _check_gene_symbol(approved_genes_same_ucsc[0], gene_symbol, db_id, db_chrom):
                    sys.stderr.write('  Found approved gene for ' + gene_symbol + ' (as ' + kind + ') with ucsc_id ' + db_id + '\n')
                    return approved_genes_same_ucsc[0].name

            # Ok, no genes with same ucsc id, or not the same chromosome for them.

            approved_genes_same_chrom = [g for g in approved_genes if g.chrom == db_chrom]

            if len(approved_genes_same_chrom) > 1:
                sys.stderr.write('  Error: multiple approved gene names for ' + gene_symbol + ' (as ' + kind + ') with chrom ' + db_chrom + ', '.join(g.name for g in approved_genes_same_ucsc) + '\n')
                return 'AMBIGUOUS'

            if len(approved_genes_same_chrom) == 1:
                g = approved_genes_same_chrom[0]
                sys.stderr.write('  Only ' + g.name + ' for ' + gene_symbol + ' (as ' + kind + ') has the same chrom ' + db_chrom + ', picking it\n')
                if _check_gene_symbol(g, gene_symbol, db_id, db_chrom):
                    return g.name
                else:
                    return 'NOT FOUND'

            if len(approved_genes_same_chrom) == 0:
                sys.stderr.write('  Error: no approved gene names for ' + gene_symbol + ' (as ' + kind + ') with same chrom ' + db_chrom + '\n')
                return 'NOT FOUND'

        if len(approved_genes) == 1:
            if _check_gene_symbol(approved_genes[0], gene_symbol, db_id, db_chrom):
                sys.stderr.write('  Found approved gene for ' + gene_symbol + ' (as ' + kind + ')\n')
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
            return None, res
        else:
            sys.stderr.write('  Finally found approved gene for ' + gene_symbol + ' (as synonym): ' + res + '\n')
            return res, None
    else:
        sys.stderr.write('  Finally found approved gene for ' + gene_symbol + ' (as prev): ' + res + '\n')
        return res, None


def _proc_ucsc(inp, out, approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym):
    not_approved_gene_names = list()

    for l in inp:
        if l and not l.startswith('#'):
            ucsc_id, ucsc_chrom, strand, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, geneSymbol = l.replace('\n', '').split('\t')

            approved_gene_symbol, status = get_approved_gene_symbol(
                approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym,
                geneSymbol, ucsc_id, ucsc_chrom)

            if approved_gene_symbol:
                txStart = exonStarts.split(',')[0]
                txEnd = exonEnds.split(',')[int(exonCount) - 1]
                out.write('\t'.join([ucsc_chrom, min(txStart, cdsStart), max(txEnd, cdsEnd),
                                     approved_gene_symbol, '.', strand, 'Gene', '.']) + '\n')
                for j, s, e in zip(range(int(exonCount)),
                   [s for s in exonStarts.split(',') if s], [
                    e for e in exonEnds.split(',') if e]):
                    if e <= cdsStart or s > cdsEnd:
                        out.write('\t'.join([ucsc_chrom, s, e, approved_gene_symbol, '.', strand, 'Exon', '.']) + '\n')
                    else:
                        if max(s, cdsStart) != s:
                            out.write('\t'.join([ucsc_chrom, s, cdsStart, approved_gene_symbol, '.', strand, 'Exon', '.']) + '\n')
                            s = cdsStart
                        if min(e, cdsEnd) != e:
                            out.write('\t'.join([ucsc_chrom, cdsEnd, e, approved_gene_symbol, '.', strand, 'Exon', '.']) + '\n')
                            e = cdsEnd
                        if s < e:  # s == e  means region of size 0
                            out.write('\t'.join([ucsc_chrom, s, e, approved_gene_symbol, '.', strand, 'CDS', '.']) + '\n')
            else:
                not_approved_gene_names.append(geneSymbol + '\t' + status)

    return not_approved_gene_names


class Gene:
    def __init__(self, name, chrom, start, end, strand, biotype, db_id, source):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.biotype = biotype
        self.db_id = db_id
        self.feature = 'Gene'
        self.source = source

        self.approved_gname = None

        self.exons = []

    def __str__(self):
        fs = [self.chrom,
              '{}'.format(self.start) if self.start else '.',
              '{}'.format(self.end) if self.end else '.',
              self.name or '.', '.', self.strand or '.',
              self.feature or '.', self.biotype or '.']
        return '\t'.join(fs) + '\n'

    def __repr__(self):
        return '{self.name} {self.chrom}:{self.start}-{self.end} {self.biotype} {self.db_id} {self.source}'.format(self=self)


class Exon:
    def __init__(self, gene, start, end, biotype, feature):
        self.gene = gene
        self.start = start
        self.end = end
        self.biotype = biotype
        self.feature = feature

    def __str__(self):
        fs = [self.gene.chrom,
              '{}'.format(self.start) if self.start else '.',
              '{}'.format(self.end) if self.end else '.',
              self.gene.name or '.', '.', self.gene.strand or '.',
              self.feature or '.', self.biotype or '.']
        return '\t'.join(fs) + '\n'


def _rm_quotes(l):
    return l[1:-1]


def is_approved_symbol(gname, approved_gene_by_name):
    if gname not in approved_gene_by_name:
        # gname2 = gname.split('.')[0]
        # if gname != gname2:
        #     if gname2 not in approved_gene_by_name:
        return False
    return True


def _proc_ensembl(inp, out, approved_gene_by_name):
    sys.stderr.write('Parsing Ensembl input...\n')
    total_lines = 0
    total_non_coding_genes = 0

    for l in inp:
        if l and not l.startswith('#'):
            chrom, _, feature, start, end, _, strand, _, props_line = l.replace('\n', '').split('\t')

            # if is_local():
            #     if chrom != '21':
            #         continue

            # total_lines += 1
            # if total_lines % 1000 == 0:
            #     sys.stderr.write(str(total_lines / 1000) + 'k lines, ' + str(len(gene_by_name)) + ' genes found\n')
            #     sys.stderr.flush()

            try:
                _prop_dict = dict((t.strip().split(' ')[0], ' '.join(t.strip().split(' ')[1:])) for t in props_line.split(';') if t.strip())
            except ValueError:
                sys.stderr.write(format_exc())
                sys.stderr.write(l)
                continue

            gene_symbol = _rm_quotes(_prop_dict['gene_name'])
            # gene_id = _rm_quotes(_prop_dict['gene_id'])
            # gene_biotype = _rm_quotes(_prop_dict['gene_biotype'])
            source = _rm_quotes(_prop_dict['gene_source'])

            # if gene_symbol == 'PTENP1':
            #     sys.stderr.write('PTENP1\n')

            # if not ALL_EXONS and gene_biotype not in [
            #     'protein_coding',
            #     'nonsense_mediated_decay',
            #     'non_stop_decay',
            #     'processed_transcript',
            #     'polymorphic_pseudogene',
            #     'sense_intronic',
            #     'sense_overlapping',
            #     'antisense',
            # ] and not any(b in gene_biotype for b in ['RNA', 'IG_', 'TR_']):
            #     total_non_coding_genes += 1
            #     continue

            if feature == 'gene' and gene_symbol in approved_gene_by_name:
                # assert gene_biotype == biotype, 'Gene: gene_biotype "' + gene_biotype + '" do not match biotype "' + biotype + '" for ' + gene_symbol
                gene = approved_gene_by_name[gene_symbol]

                start, end = int(start) - 1, int(end)
                if int(end) <= int(start):
                    sys.stderr.write('Error: start > end: ' + l + '\n')
                    continue
                length = end - start

                if gene.source is None:
                    gene.length = length

                else:
                    prev_source = gene.source

                    if source != prev_source:
                        sys.stderr.write('    Duplicated gene in different databases:\n')
                        sys.stderr.write('        This: ' + source + ', ' + str(length) + '\n')
                        sys.stderr.write('        Prev: ' + prev_source + ', ' + str(gene.length) + '\n')
                        # answer = raw_input('Which one to pick? This (1), prev (2), longest (Enter): ')
                        #
                        # if answer == '1' or answer == '' and gene.end - gene.start > prev_gene.end - prev_gene.start:
                        #     del gene_by_name[prev_gene.name]
                        #     del gene_by_id[prev_gene.db_id]
                        #
                        # else:
                        #     continue

                        if source == 'ensembl' or prev_source == 'havana':
                            gene.source = source
                            gene.length = length
                            sys.stderr.write('        Picking up this one.\n')

                        if prev_source == 'ensembl' or source == 'havana':
                            sys.stderr.write('        Picking up previous one.\n')
                            continue

                    else:
                        sys.stderr.write('    MultiGene gene in ' + gene.source + ':\n')
                        sys.stderr.write('        ' + gene.__repr__() + '\n')
                        gene.length += length
                        continue

                    sys.stderr.write('\n')

    sys.stderr.write('\n')
    sys.stderr.write(
        'Processed ' +
        str(total_lines) + ' lines, ' +
        str(total_non_coding_genes) + ' non-coding genes skipped, ' +
        str(len(approved_gene_by_name)) + ' coding genes found\n')
    sys.stderr.write('\n')

    not_found_genes_total = 0
    found_genes_total = 0

    for g in approved_gene_by_name.values():
        if g.length is None:
            not_found_genes_total += 1
            sys.stderr.write('  gene ' + g.name + ' not found in Ensemble\n')
        else:
            found_genes_total += 1
            out.write(str(g) + '\n')
        # for e in g.exons:
        #     f.write('\t' + str(e) + '\n')

    sys.stderr.write('\n')
    sys.stderr.write('Found genes in Ensembl ' + str(found_genes_total) + '\n')
    sys.stderr.write('Not-found genes in Ensembl ' + str(not_found_genes_total) + '\n')
    return

    not_approved_gene_names = dict()
    gene_after_approving_by_name = OrderedDict()
    total_approved = 0
    total_not_approved = 0
    j = 0
    for g in approved_gene_by_name.values():
        if len(g.exons) == 0:
            continue

        if not DO_APPROVE:
            gene_after_approving_by_name[g.name] = g
        if not approved_gene_by_name or is_approved_symbol(g.name, approved_gene_by_name):
            if DO_APPROVE:
                gene_after_approving_by_name[g.name] = g
            total_approved += 1
        else:
            total_not_approved += 1

        j += 1
        if j % 1000 == 0:
            sys.stderr.write('processed ' + str(j / 1000) + 'k genes...\n')

    sys.stderr.write('-----\n')
    sys.stderr.write('Total: ' + str(j) + '\n')
    if approved_gene_by_name:
        sys.stderr.write('Total approved: ' + str(total_approved) + '\n')
        sys.stderr.write('Total not approved: ' + str(total_not_approved) + '\n')
    sys.stderr.write('\n')
    sys.stderr.write('Saving genes...\n')

    gene_features = 0
    features_counter = defaultdict(int)
    biotypes_counter = defaultdict(int)
    no_exon_gene_num = 0

    for g in gene_after_approving_by_name.values():
        if len(g.exons) == 0:
            no_exon_gene_num += 1
        else:
            out.write(str(g))

            gene_features += 1
            features_counter[g.feature] += 1
            biotypes_counter[g.biotype] += 1

            for e in g.exons:
                features_counter[e.feature] += 1

                if e.feature == 'exon': e.feature = 'Exon'
                elif e.feature == 'stop_codon': e.feature = 'CDS'
                else: e.feature = e.feature[0].upper() + e.feature[1:]

                out.write(str(e))

    sys.stderr.write('Skipped {} genes with no sub-features.\n'.format(no_exon_gene_num))
    sys.stderr.write('Saved {} genes, including:\n'.format(gene_features))
    sys.stderr.write('    Gene: {}\n'.format(features_counter['Gene']))
    sys.stderr.write('    Multi_Gene: {}\n'.format(features_counter['Multi_Gene']))
    sys.stderr.write('\n')

    sys.stderr.write('Out of total: {} protein coding genes, {} ncRNA genes, including:\n'.format(
        biotypes_counter['protein_coding'], sum(biotypes_counter.values()) - biotypes_counter['protein_coding']))
    for bt, cnt in biotypes_counter.items():
        if bt != 'protein_coding':
            sys.stderr.write('    ' + bt + ': ' + str(cnt) + '\n')

    sys.stderr.write('\n')
    if ALL_EXONS:
        sys.stderr.write('Found {} exons.\n'.format(features_counter['exon']))
    else:
        sys.stderr.write('Also found {} CDS, {} stop codons, and {} ncRNA exons.\n'.format(
            features_counter['CDS'], features_counter['stop_codon'], features_counter['exon']))

    return not_approved_gene_names


if __name__ == '__main__':
    main()
