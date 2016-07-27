#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

import sys
import re
from collections import defaultdict
from os.path import join, basename, isdir

import source
from source import info
from source.clinical_reporting.utils import get_key_or_target_bed_genes, SVEvent
from source.file_utils import verify_file, add_suffix, file_transaction, adjust_system_path
from source.bcbio.bcbio_structure import bcbio_summary_script_proc_params, BCBioStructure
from source.logger import critical, warn, err, step_greetings
from source.prepare_args_and_cnf import check_genome_resources
from source.clinical_reporting.known_sv import fusions as known_fusions
from source.targetcov.Region import get_chrom_order
from source.utils import is_uk
from source.utils import is_us
from source.webserver import exposing
from source.webserver.exposing import symlink_to_ngs, local_symlink


class OncoprintMutation:
    def __init__(self, chrom, pos, gene):
        self.gene = gene
        self.chrom = chrom
        self.pos = pos
        self.type = None
        self.aa_change = None
        self.cdna_change = None
        self.status = None  # germline/somatic
        self.signif = None
        self.zygosity = None

        self.depth = None
        self.freq = None

    def __str__(self):
        return str(self.gene) + ' ' + str(self.chrom) + ':' + \
               str(self.pos) + ' ' + str(self.aa_change) + ' ' + str(self.signif)

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash((self.chrom, self.pos, self.gene, self.aa_change))

    def get_key(self):
        return self.chrom, self.pos, self.gene, self.aa_change


class OncoprintSeq2CEvent():
    def __init__(self, gene=None, exons=None, copy_number=None, ratio=None, cnv_type=None):
        self.gene = gene
        self.copy_number = copy_number
        self.exons = exons
        self.ratio = ratio
        self.cnv_type = cnv_type # loss, amplification


def proc_args():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = bcbio_summary_script_proc_params(
        BCBioStructure.seq2c_name,
        extra_opts=[
           (['--bed', '--capture', '--amplicons'], dict(
                dest='bed'
           ))
        ],
    )
    return cnf, bcbio_structure


def main():
    cnf, bcbio_structure = proc_args()
    create_oncoprints_link(cnf, bcbio_structure)


def create_oncoprints_link(cnf, bcbio_structure, project_name=None):
    if is_us(): loc = exposing.us
    # elif is_uk(): loc = exposing.uk
    else:
        loc = exposing.local
        return None

    step_greetings('Creating Oncoprints link')
    zhongwu_data_query_dirpath = '/home/kdld047/public_html/cgi-bin/TS'
    if not isdir(zhongwu_data_query_dirpath):
        warn('Data Query directory ' + zhongwu_data_query_dirpath + ' does not exists.')
        return None

    clinical_report_caller = \
        bcbio_structure.variant_callers.get('vardict') or \
        bcbio_structure.variant_callers.get('vardict-java')
    if not clinical_report_caller:
        critical('No vardict or vardict-java variant caller in ' + str(bcbio_structure.variant_callers.keys()))
    vardict_txt_fname = source.mut_fname_template.format(caller_name=clinical_report_caller.name)
    vardict_txt_fpath = join(bcbio_structure.var_dirpath, vardict_txt_fname)
    cnf.mutations_fpath = add_suffix(vardict_txt_fpath, source.mut_pass_suffix)

    cnf.seq2c_tsv_fpath = bcbio_structure.seq2c_fpath

    samples = sorted(bcbio_structure.samples)
    cnf.project_name = project_name or bcbio_structure.project_name or basename(cnf.output_dir)
    study_name = re.sub('[\.\-:&]', '_', cnf.project_name)

    check_genome_resources(cnf)

    data_query_dirpath = join(loc.dirpath, 'DataQueryTool')

    data_fpath = join(zhongwu_data_query_dirpath, study_name + '.data.txt')
    info_fpath = join(zhongwu_data_query_dirpath, study_name + '.info.txt')
    altered_genes = print_data_txt(cnf, cnf.mutations_fpath, cnf.seq2c_tsv_fpath, samples, data_fpath)
    if not altered_genes:
        err('No altered genes!')
        return None

    print_info_txt(cnf, samples, info_fpath)

    data_ext_fpath = data_fpath.replace('/home/', '/users/')
    info_ext_fpath = info_fpath.replace('/home/', '/users/')

    # optional:
    data_symlink = join(data_query_dirpath, study_name + '.data.txt')
    info_symlink = join(data_query_dirpath, study_name + '.info.txt')
    (symlink_to_ngs if is_us() else local_symlink)(data_ext_fpath, data_symlink)
    (symlink_to_ngs if is_us() else local_symlink)(info_ext_fpath, info_symlink)

    properties_fpath = join(zhongwu_data_query_dirpath, 'DataQuery.properties')
    add_data_query_properties(cnf, study_name, properties_fpath, data_ext_fpath, info_ext_fpath)

    genes = '%0D%0A'.join(altered_genes)
    data_query_url = join(loc.website_url_base, 'DataQueryTool', 'DataQuery.pl?'
        'analysis=oncoprint&'
        'study={study_name}&'
        'gene={genes}&'
        'order=on&'
        'freq=50&'
        'nocheckgenes=true&'
        'submit=Submit'
        .format(**locals()))

    info()
    info('Information about study was added in Data Query Tool, URL is ' + data_query_url)
    return data_query_url


def print_data_txt(cnf, mutations_fpath, seq2c_tsv_fpath, samples, data_fpath):
    bed_fpath = verify_file(cnf.bed, is_critical=False) if cnf.bed else None
    key_gene_by_name_chrom, _ = get_key_or_target_bed_genes(bed_fpath, verify_file(adjust_system_path(cnf.key_genes), 'key genes'))
    key_genes = [g for (g, c) in key_gene_by_name_chrom]

    altered_genes = set()

    mut_by_samples, altered_genes = parse_mutations(mutations_fpath, altered_genes, key_genes)
    seq2c_events_by_sample, altered_genes = parse_seq2c(seq2c_tsv_fpath, altered_genes, key_genes)
    sv_events_by_samples, altered_genes = parse_sv_files(cnf, samples, altered_genes, key_genes)

    with file_transaction(cnf.work_dir, data_fpath) as tx:
        with open(tx, 'w') as out_f:
            out_f.write('SAMPLE ID\tANALYSIS FILE LOCATION\tVARIANT-TYPE\tGENE\tSOMATIC STATUS/FUNCTIONAL IMPACT\tSV-PROTEIN-CHANGE\t'
                        'SV-CDS-CHANGE\tSV-GENOME-POSITION\tSV-COVERAGE\tSV-PERCENT-READS\tCNA-COPY-NUMBER\tCNA-EXONS\tCNA-RATIO\t'
                        'CNA-TYPE\tREARR-GENE1\tREARR-GENE2\tREARR-DESCRIPTION\tREARR-IN-FRAME?\tREARR-POS1\tREARR-POS2\tREARR-NUMBER-OF-READS\n')
            for sample, muts in mut_by_samples.iteritems():
                for mut in muts:
                    out_f.write('\t'.join([sample, '', 'short-variant', mut.gene, mut.signif, mut.aa_change, mut.cdna_change,
                                           mut.chrom + ':' + mut.pos, str(mut.depth), str(mut.freq * 100),
                                           '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', mut.type]) + '\n')
            for sample, events in seq2c_events_by_sample.iteritems():
                for event in events:
                    out_f.write('\t'.join([sample, '', 'copy-number-alteration', event.gene, 'NA', '-', '-',
                                           '-', '-', '-', str(event.copy_number), event.exons, str(event.ratio), event.cnv_type,
                                           '-', '-', '-', '-', '-', '-', '-', 'Deletion' if event.cnv_type == 'loss' else 'Amplification']) + '\n')
            for sample, events in sv_events_by_samples.iteritems():
                for event in events:
                    for ann in event.key_annotations:
                        count_reads = sum(int(r) for r in event.split_read_support) + sum(int(r) for r in event.paired_end_support)
                        out_f.write('\t'.join([sample, '', 'rearrangement', event.known_gene, 'known' if ann.known else 'unknown',
                                               '-', '-', '-', '-', '-', '-', '-', '-', '-', ann.genes[0], ann.genes[1],
                                               'fusion', event.chrom + ':' + str(event.start), str(event.end) if event.end else '',
                                               '-', str(count_reads), 'Rearrangement']) + '\n')

    return altered_genes


def parse_mutations(mutations_fpath, altered_genes, key_genes):
    mut_by_samples = defaultdict(list)
    
    if not mutations_fpath or not verify_file(mutations_fpath):
        return mut_by_samples, altered_genes

    info('Parsing mutations from ' + mutations_fpath)
    
    sample_col = None
    chr_col = None
    pos_col = None
    type_col = None
    allele_freq_col = None
    gene_col = None
    depth_col = None
    aa_chg_col = None
    cdna_chg_col = None
    status_col = None
    signif_col = None
    incidentalome_col = None

    stop_gain_pattern = re.compile('^[A-Z]+\d+\*')
    fs_pattern = re.compile('^[A-Z]+(\d+)fs')
    aa_chg_pattern = re.compile('^([A-Z]\d+)[A-Z]$')
    
    with open(mutations_fpath) as txt:
        for i, l in enumerate(txt):
            l = l.replace('\n', '')
            if not l:
                continue
            if i == 0:
                header = l.split('\t')
                sample_col = header.index('Sample')
                chr_col = header.index('Chr')
                pos_col = header.index('Start')
                type_col = header.index('Type')
                allele_freq_col = header.index('AlleleFreq')
                gene_col = header.index('Gene')
                aa_chg_col = header.index('Amino_Acid_Change')
                cdna_chg_col = header.index('cDNA_Change')
                depth_col = header.index('Depth')
                if 'Status' in header:
                    status_col = header.index('Status')
                if 'Significance' in header:
                    signif_col = header.index('Significance')
                else:
                    signif_col = len(header) - header[::-1].index('Status') - 1  # get last index of status
                if 'Incidentalome' in header:
                    incidentalome_col = header.index('Incidentalome')
                continue
            fs = l.replace('\n', '').split('\t')
            sample, gene, chrom, pos, type_ = fs[sample_col], fs[gene_col], fs[chr_col], fs[pos_col], fs[type_col]
            if gene not in key_genes:
                continue
            mut = OncoprintMutation(chrom, pos, gene)
            mut.aa_change, mut.cdna_change, mut.depth, mut.freq = fs[aa_chg_col], fs[cdna_chg_col], fs[depth_col], float(fs[allele_freq_col])
            mut.status = fs[status_col] if status_col is not None else None
            mut.signif = fs[signif_col] if signif_col is not None else None
            incidentalome_reason = fs[incidentalome_col] if incidentalome_col is not None else None
            if incidentalome_reason:
                continue
            mut.type = 'Known' if mut.signif != 'unknown' else 'Unknown'
            if 'splice' in type_:
                mut.type = 'Splice'
            elif stop_gain_pattern.match(mut.aa_change):
                mut.type = 'Trunc/FS'
            elif fs_pattern.match(mut.aa_change):
                mut.type = 'Trunc/FS'
            elif mut.aa_change.startswith('-'):
                mut.type = 'Trunc/FS'
            elif 'missense' in type_:
                mut.type += '-Missense'
            elif 'ins' in mut.aa_change or 'del' in mut.aa_change:
                mut.type += '-Indel'
            else:
                mut.type += '-Other'
            mut_by_samples[sample].append(mut)
            altered_genes.add(gene)

    return mut_by_samples, altered_genes


def parse_seq2c(seq2c_tsv_fpath, altered_genes, key_genes):
    seq2c_events_by_sample = defaultdict(list)
    
    if not seq2c_tsv_fpath or not verify_file(seq2c_tsv_fpath):
        return seq2c_events_by_sample, altered_genes

    info('Parsing Seq2C from ' + seq2c_tsv_fpath)

    if seq2c_tsv_fpath and verify_file(seq2c_tsv_fpath):
        with open(seq2c_tsv_fpath) as f_inp:
            for i, l in enumerate(f_inp):
                if i == 0: continue
                fs = l.replace('\n', '').split('\t')
                sname, gname = fs[0], fs[1]
                if gname not in key_genes: continue

                sname, gname, chrom, start, end, length, log2r, sig, fragment, amp_del, ab_seg, total_seg, \
                    ab_log2r, log2r_diff, ab_seg_loc, ab_samples, ab_samples_pcnt = fs[:17]
                if not amp_del:
                    continue
                if gname not in key_genes:
                    continue
                if fragment == 'BP':
                    exons = str(ab_seg) + ' of ' + total_seg
                    copy_number = round(2 ** float(ab_log2r) * 2, 2)
                    if copy_number > 2:
                        copy_number = round(copy_number, 1)
                else:
                    exons = 'Whole'
                    if amp_del == 'Amp':
                        ab_log2r = copy_number = 'AMP'
                    else:
                        ab_log2r = copy_number = 'HOMDEL'

                cnv_type = 'amplification' if amp_del == 'Amp' else 'loss'
                event = OncoprintSeq2CEvent(
                    gene=gname,
                    copy_number=str(copy_number),
                    exons=exons,
                    ratio=ab_log2r,
                    cnv_type=cnv_type)
                seq2c_events_by_sample[sname].append(event)
                altered_genes.add(gname)
                
    return seq2c_events_by_sample, altered_genes


def parse_sv_files(cnf, samples, altered_genes, key_genes):
    sv_events_by_samples = defaultdict(set)
    sv_fpaths = [sample.find_sv_fpath() for sample in samples]
    sv_fpaths = [f for f in sv_fpaths if f]

    if not sv_fpaths:
        return sv_events_by_samples, altered_genes
    
    sorted_known_fusions = [sorted(p) for p in known_fusions['homo_sapiens']]

    chr_order = get_chrom_order(cnf)

    for sv_fpath in sv_fpaths:
        info('Parsing prioritized SV from ' + sv_fpath)
        sample_col = None
        known_col = None
        with open(sv_fpath) as f:
            header_rows = []
            for i, l in enumerate(f):
                fs = l.strip().split('\t')
                if i == 0:
                    header_rows = fs  # caller  sample  chrom  start  end  svtype  known  end_gene  lof  annotation  split_read_support  paired_end_support
                    sample_col = header_rows.index('sample')
                    known_col = header_rows.index('known')
                else:
                    event = SVEvent.parse_sv_event(chr_order, **dict(zip(header_rows, fs)))
                    sample = fs[sample_col]
                    if event:
                        for annotation in event.annotations:
                            if event.is_fusion():
                                if sorted(annotation.genes) in sorted_known_fusions:
                                    annotation.known = True
                                key_altered_genes = [g for g in annotation.genes if g in key_genes]
                                if annotation.effect == 'FUSION' and key_altered_genes:
                                    event.key_annotations.add(annotation)
                                    event.supplementary = '-with-' in fs[known_col]
                                    sv_events_by_samples[sample].add(event)
                                    for g in key_altered_genes:
                                        altered_genes.add(g)

    for sample, events in sv_events_by_samples.iteritems():  # combine two annotations of fusion in one
        suppl_events = {e.mate_id: e for e in events if e.supplementary}
        main_events = [e for e in events if not e.supplementary]
        for event in main_events:
            if event.id not in suppl_events:
                continue
            suppl_event = suppl_events[event.id]
            event.end = suppl_event.chrom + ':' + str(suppl_event.start)
            for ann in event.key_annotations:
                if ann.genes[0] != event.known_gene:
                    ann.genes[0], ann.genes[1] = ann.genes[1], ann.genes[0]
        sv_events_by_samples[sample] = main_events

    return sv_events_by_samples, altered_genes


def print_info_txt(cnf, samples, info_fpath):
    with file_transaction(cnf.work_dir, info_fpath) as tx:
        with open(tx, 'w') as out_f:
            out_f.write('Sample\n')
            out_f.write('Type\n')
            for sample in samples:
                out_f.write(sample.name + '\n')


# DATA_QUERY_LINK = 'http://ngs.usbod.astrazeneca.net/DataQueryTool/DataQuery.pl'

def add_data_query_properties(cnf, study_name, properties_fpath, data_fpath, info_fpath):
    # modify properties
    properties_lines = []
    text_to_add = None

    lines = open(properties_fpath).read().split('\n')
    for l in lines:
        l = l.strip()
        if 'studies=' in l:
            studies = l.split('=')[1].split(';')
            if study_name not in studies:
                l += ';' + study_name
        if 'study2desc' in l:
            text_to_add = '{study_name} => "{study_name}", \\'.format(**locals())
        if 'study2data' in l:
            text_to_add = '{study_name} => "{data_fpath}", \\'.format(**locals())
        if 'study2info' in l:
            text_to_add = '{study_name} => "{info_fpath}", \\'.format(**locals())
        if study_name + ' => ' in l:
            info(l.strip() + ' already present in properties, removing it.')
            continue
        if l == '}' and text_to_add and text_to_add not in properties_lines:
            properties_lines.append(text_to_add)
            text_to_add = None
        properties_lines.append(l)

    with file_transaction(cnf.work_dir, properties_fpath) as tx:
        with open(tx, 'w') as out:
            for l in properties_lines:
                out.write(l + '\n')


if __name__ == '__main__':
    main()
