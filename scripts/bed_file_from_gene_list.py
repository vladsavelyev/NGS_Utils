#!/usr/bin/env python
import sys
from collections import defaultdict

from Utils.tools import read_approved_genes, get_approved_gene_symbol


def main():
    if len(sys.argv) < 2:
        log('Usage: ' + __file__ + ' genes.txt [regions.bed] [HGNC_gene_synonyms.txt] [--not-check] > regions_for_genes.txt')
        log('    the script filters regions to have gene names in genes.txt')
        sys.exit(1)

    genes_list_fpath = sys.argv[1]
    exons_fpath = '/ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/RefSeq/RefSeq_CDS.hg19.bed'
    synonyms_fpath = '/ngs/reference_data/genomes/Hsapiens/common/HGNC_gene_synonyms.txt'
    if len(sys.argv) > 2: exons_fpath = sys.argv[2]
    if len(sys.argv) > 3: synonyms_fpath = sys.argv[3]
    not_check = len(sys.argv) > 4

    approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym = \
        read_approved_genes(synonyms_fpath)

    gnames = set()
    duplicated_genes = list()
    with open(genes_list_fpath) as f:
        for i, l in enumerate(f):
            gname = l.strip()
            sys.stderr.write(str(i + 1) + ' ' + gname + ' ')

            if gname in gnames:
                duplicated_genes.append(gname)
                sys.stderr.write(' duplicated gene\n')
                continue

            gnames.add(gname)
            sys.stderr.write('\n')

    log('Total unique gene names: ' + str(len(gnames)))
    if duplicated_genes:
        log(str(len(duplicated_genes)) + ' duplicated: ' + ', '.join(duplicated_genes))
    log()

    genes_found_in_ref, lines = exons_for_gene_list(exons_fpath, gnames)
    # log(str(len(doublicated_after_approving)) + ' duplicated after approval: ' + ', '.join(doublicated_after_approving))
    # log('Not found in HGNC ' + str(len(not_found_in_HGNC_genes)) + ' genes' +
    #     (': ' + ', '.join(not_found_in_HGNC_genes) if not_found_in_HGNC_genes else ''))
    log('Found in Exons DB ' + str(len(genes_found_in_ref)) + ' genes out of ' + str(len(gnames)))
    not_found_in_reference = gnames - genes_found_in_ref
    log('Not found in Exons DB ' + str(len(not_found_in_reference)) + ': ' + ', '.join(not_found_in_reference))

    doublicated_after_approving_gnames = list()
    doublicated_after_approving_approved = list()
    original_by_approved = defaultdict(list)
    corrected_genes = set()
    corrected_not_found_in_HGNC_genes = set()
    not_corrected = set()
    if not_found_in_reference:
        log()
        log('Correcting genes not found in Exons DB')
        for gname in not_found_in_reference:
            log(gname)

            approved_gname, status = get_approved_gene_symbol(
                approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym, gname,
                indent='    ')

            if not approved_gname:
                if '.' in gname:
                    gname2 = gname.split('.')[0]
                    log('    Not found with dot, trying without dot, as ' + gname2 + '...')
                    approved_gname, status2 = get_approved_gene_symbol(
                        approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym, gname2,
                        indent='    ')

            if approved_gname:
                if gname == approved_gname:
                    log('    Gene name ' + gname + ' is already approved. Skipping.')
                    not_corrected.add(gname)
                else:
                    log('    Found as ' + approved_gname)

                    if approved_gname in gnames:
                        log('    Approved version for ' + gname + ' is already met as ' + approved_gname +
                                         ' in gnames')
                        original_by_approved[approved_gname].append(gname)
                        doublicated_after_approving_gnames.append(gname + '->' + approved_gname)
                        continue

                    if approved_gname in original_by_approved:
                        log('    Approved version for ' + gname + ' is already met as ' + approved_gname +
                                         ' for ' + ', '.join(original_by_approved[approved_gname]))
                        original_by_approved[approved_gname].append(gname)
                        doublicated_after_approving_approved.append(gname + '->' + approved_gname)
                        continue

                    log('    Found unique corrected approved gene name ' + approved_gname)
                    original_by_approved[approved_gname].append(gname)
                    corrected_genes.add(approved_gname)

            else:
                log('    ' + gname + ' not found in HGNC.')
                corrected_not_found_in_HGNC_genes.add(gname)
            log()

    log('-' * 70)
    log('Total unique gene names: ' + str(len(gnames)))
    if duplicated_genes:
        log('    ' + str(len(duplicated_genes)) + ' duplicated: ' + ', '.join(duplicated_genes))

    log()
    log('Searching exons for the gene list: first iteration')
    log('    Found ' + str(len(lines)) + ' exons for ' + str(len(genes_found_in_ref)) + ' genes (out of ' + str(len(gnames)) + ')')
    not_found_in_reference = gnames - genes_found_in_ref
    log('    Not found in Exons DB ' + str(len(not_found_in_reference)) + ': ' + ', '.join(not_found_in_reference))

    corrected_genes_found_in_ref = set()
    corrected_not_found_in_reference = set()
    if not_found_in_reference:
        log()
        log('Checking failed gene names with HGNC')
        log('    ' + str(len(doublicated_after_approving_gnames)) + ' approved names are already accounted for in the first iteration: ' + ', '.join(doublicated_after_approving_gnames))
        log('    ' + str(len(doublicated_after_approving_approved)) + ' duplicated approved names' + ((': ' + ', '.join(doublicated_after_approving_approved)) if doublicated_after_approving_approved else ''))
        log('    ' + str(len(corrected_not_found_in_HGNC_genes)) + ' not found in HGNC' + (': ' + ', '.join(corrected_not_found_in_HGNC_genes) if corrected_not_found_in_HGNC_genes else ''))
        log('    ' + str(sum([len(original_by_approved[c]) for c in corrected_genes])) + ' gene names corrected by HGNC into new names: ' + ', '.join([':'.join(original_by_approved[c]) + '->' + c for c in corrected_genes]))
        log('    ' + str(len(not_corrected)) + ' not corrected gene names.')

        log('Writing exons for the approved version of gene names originally not found in the Exons DB')
        corrected_genes_found_in_ref, lines_2 = exons_for_gene_list(exons_fpath, corrected_genes)
        corrected_not_found_in_reference = corrected_genes - corrected_genes_found_in_ref
        log('    Found ' + str(len(lines_2)) + ' new exons for ' + str(len(corrected_genes_found_in_ref)) + ' genes (out of ' + str(len(corrected_genes)) + ')')
        log('    Out of corrected, not found in Exons DB ' + str(len(corrected_not_found_in_reference)) + ((': ' + ', '.join(corrected_not_found_in_reference)) if corrected_not_found_in_reference else ''))

    log('')
    total_uniq_genes_found = len(gnames) - len(corrected_not_found_in_HGNC_genes) - len(corrected_not_found_in_reference)
    log('Finally writing exons for ' + str(total_uniq_genes_found) + ' unique gene name records, ' +
        str(len(genes_found_in_ref | corrected_genes_found_in_ref)) + ' unique genes after correcting with HGNC')
    genes_found_in_ref, lines = exons_for_gene_list(exons_fpath, genes_found_in_ref | corrected_genes_found_in_ref)
    log('    Written ' + str(len(lines)) + ' lines for ' + str(len(genes_found_in_ref)) + ' genes')
    missed_genes = not_corrected | corrected_not_found_in_HGNC_genes
    log('    Missed ' + str(len(missed_genes)) + ' genes: ' + ', '.join(missed_genes))

    for l in lines:
        print(l)


def exons_for_gene_list(exons_fpath, gnames):
    genes_found_in_ref = set()
    lines = []

    with open(exons_fpath) as f:
        for l in f:
            if not l.startswith('#'):
                fs = l.strip().split('\t')
                if len(fs) < 7:
                    pass
                else:
                    feature = fs[6]
                    # if feature in ['Gene', 'Multi_Gene']:
                    #     pass
                    # else:
                    gname = fs[3]
                    if gname not in gnames:
                        pass
                        # genes_not_found_in_ref.add(gname)
                    else:
                        genes_found_in_ref.add(gname)
                        lines.append('\t'.join(fs))

    return genes_found_in_ref, lines


def log(msg=''):
    sys.stderr.write(msg + '\n')


if __name__ == '__main__':
    main()
