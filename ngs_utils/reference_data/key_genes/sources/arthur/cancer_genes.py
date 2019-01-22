#!/usr/bin/env python

from collections import defaultdict
from datetime import date
import sys

def get_genes(gene_file):
    genes = [ line.rstrip().split('\t')[0] for line in open(gene_file) ]
    genes = set(genes)

    return genes


clinical_panels = ["PMCC", "TS500", "MSKCC_IMPACT.curated", "MSKCC_HEME", "FoundationOne", "FoundationHeme", "TEMPUS"]
database_genes = ["CancerGeneCensus_Tier1", "OncoKB-Annotated"]
research_genes = ["CancerGeneCensus_Tier2", "az300", "familial_cancer", "TSGene2-tumour-suppressors", "TSGene2-oncogenes"]

gene_set_order = clinical_panels + database_genes + research_genes

gene_sets = {}
all_genes = set([])
for gene_set_name in gene_set_order:
    gene_set_file = gene_set_name + '.genes'
    gene_sets[gene_set_name] = get_genes(gene_set_file)
    all_genes |= gene_sets[gene_set_name]


total = 0
ensemble_gene_list = 0
in_any_database = 0
in_any_panel = 0
unique_gene_in_panel = dict( [ (p, 0) for p in gene_set_order ] )

outfile1 = 'cancer_genes.{:04d}{:02d}{:02d}.txt'.format(*(date.today().timetuple()[:3]))
outfile2 = 'cancer_genes.{:04d}{:02d}{:02d}.genes'.format(*(date.today().timetuple()[:3]))
template = '''list(query = intersects, params = list({}), color = "{}", active = T)'''

with open(outfile1, 'w') as out1, open(outfile2, 'w') as out2:
    header = ['Gene'] + gene_set_order
    print('\t'.join(header), file=out1)

    membership_counter = defaultdict(list)
    for gene in sorted(list(all_genes)):
        panel_sources = []
        database_sources = []
        gene_list_sources = []
        membership = [0] * len(gene_set_order)

        for idx, gene_set_name in enumerate(gene_set_order):
            if gene in gene_sets[gene_set_name]:
                membership[idx] = 1

                if gene_set_name in clinical_panels:
                    panel_sources.append(gene_set_name)
                elif gene_set_name == 'CancerGeneCensus_Tier1' or gene_set_name == 'OncoKB-Annotated':
                    database_sources.append(gene_set_name)
                else:
                    gene_list_sources.append(gene_set_name)

        union_src = panel_sources + database_sources + gene_list_sources
        src = '|'.join(union_src)
        if len(panel_sources) > 1:
            in_any_panel += 1
            total += 1
            if len(union_src) == 1:
                unique_gene_in_panel[union_src[0]] += 1
            print(gene, src, sep='\t', file=out2)
        elif len(database_sources) > 0:
            in_any_database += 1
            total += 1
            if len(union_src) == 1:
                unique_gene_in_panel[union_src[0]] += 1
            print(gene, src, sep='\t', file=out2)
        else:
            if len(gene_list_sources) >= 2:
                total += 1
                ensemble_gene_list += 1
                if len(union_src) == 1:
                    unique_gene_in_panel[union_src[0]] += 1
                print(gene, src, sep='\t', file=out2)

        print(gene, '\t'.join([str(v) for v in membership]), sep='\t', file=out1)

        

    print('{} genes in cancer gene list'.format(total), file=sys.stderr)
    print('{} genes in panel\n{} genes in database\n{} genes in ensemble (2+) lists'.format(in_any_panel, in_any_database, ensemble_gene_list), file=sys.stderr)
    for p in gene_set_order:
        print('{}\t{}'.format(p, unique_gene_in_panel[p]), file=sys.stderr)


#print('In panel: {}'.format(in_any_panel), file=sys.stderr)
#print('Ensemble gene list: {}'.format(ensemble_gene_list), file=sys.stderr)
#
#
#print('UpsetR colouring:', file=sys.stderr)
#popular_membership = sorted(membership_counter.items(), key=lambda x: len(x[1]), reverse=True)
#colour_queries = []
#for scores, genes in popular_membership:
#    if len(genes) < 4:
#        break
#    sources = [gene_set_order[i] for i, p in enumerate(scores) if p > 0]
#
#    total_score = sum(scores)
#    if total_score >= 3:
#        if sum(scores[:3]) > 0:
#            colour = 'blue'
#        else:
#            colour = 'orange'
#        colour_queries.append(template.format(', '.join(['"{}"'.format(s) for s in sources]), colour))
#print('queries = list({})'.format(', \n'.join(colour_queries)), file=sys.stderr)
#
