#!/usr/bin/env python

import copy
import os
import shutil
import subprocess
import sys
import tempfile
from optparse import OptionParser
from os.path import abspath, join, splitext
from os.path import exists, basename

from source.config import Config
from source.file_utils import add_suffix, adjust_path
from source.logger import critical, info
from source.prepare_args_and_cnf import determine_sys_cnf
from source.targetcov.Region import SortableByChrom
from source.utils import get_chr_lengths_from_seq

from Utils.tools import annotate_bed

""" Input: Any BED file
    Output:
        BED file with regions from input and with exactly 4 columns (or exactly 8 columns for BEDs with primers info).
        The 4th column ("gene symbol") is generated with annotate_bed.py script. If no annotation is found
        gene symbol is set to "not_a_gene_%d". If more than one annotation is found for region only one is remained
        (priory from highest to lowest: key gene, approved gene, first gene).
        Output BED is sorted using by chromosome name -> start -> end. Run standardize_bed.py --help for details about
        options.

    Usage: python standardize_bed.py [options] Input_BED_file -o Standardized_BED_file
"""


def _read_args(args_list):
    options = [
        # (['-k', '--key-genes'], dict(
        #     dest='key_genes_fpath',
        #     help='list of key genes (they are at top priority when choosing one of multiple annotations)',
        #     default='/ngs/reference_data/genomes/Hsapiens/common/az_key_genes.300.txt')
        #  ),
        # (['-a', '--approved-genes'], dict(
        #     dest='approved_genes_fpath',
        #     help='list of HGNC approved genes (they are preferable when choosing one of multiple annotations)',
        #     default='/ngs/reference_data/genomes/Hsapiens/common/HGNC_gene_synonyms.txt')
        #  ),
        # (['-e', '--ensembl-bed'], dict(
        #     dest='ensembl_bed_fpath',
        #     help='reference BED file for annotation (Ensembl)')
        #  ),
        # (['-r', '--refseq-bed'], dict(
        #     dest='refseq_bed_fpath',
        #     help='reference BED file for annotation (RefSeq)')
        #  ),
        # (['-b', '--bedtools'], dict(
        #     dest='bedtools',
        #     help='path to bedtools',
        #     default='bedtools')
        #  ),
        (['-o', '--output-bed'], dict(
            dest='output_fpath')
         ),
        (['--debug'], dict(
            dest='debug',
            help='run in a debug more (verbose output, keeping of temporary files)',
            default=False,
            action='store_true')
         ),
        (['--output-hg'], dict(
            dest='output_hg',
            help='output chromosome names in hg-style (chrM, chr1, .., chr22, chrX, chrY)',
            default=False,
            action='store_true')
         ),
        (['--output-grch'], dict(
            dest='output_grch',
            help='output chromosome names in GRCh-style (1, .., 22, X, Y, MT)',
            default=False,
            action='store_true')
         ),
        (['-g', '--genome'], dict(
            dest='genome',
            default='hg19')
         ),
    ]

    parser = OptionParser(usage='usage: %prog [options] Input_BED_file -o Standardized_BED_file',
                          description='Scripts outputs a standardized version of input BED file. '
                                      'Standardized BED: 1) has 4 or 8 fields (for BEDs with primer info);'
                                      ' 2) has HGNC approved symbol in forth column if annotation is '
                                      'possible and not_a_gene_X otherwise;'
                                      ' 3) is sorted based on chromosome name -> start -> end;'
                                      ' 4) has no duplicated regions (regions with the same chromosome, start and end), '
                                      'the only exception is _CONTROL_ regions.')
    for args, kwargs in options:
        parser.add_option(*args, **kwargs)
    (opts, args) = parser.parse_args(args_list)

    if len(args) != 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    cnf = Config(opts.__dict__, determine_sys_cnf(opts), {})

    work_dirpath = tempfile.mkdtemp()
    info('Creating a temporary working directory ' + work_dirpath)
    if not exists(work_dirpath):
        os.mkdir(work_dirpath)

    input_bed_fpath = abspath(args[0])
    info('Input: ' + input_bed_fpath)

    output_bed_fpath = adjust_path(cnf.output_fpath)
    info('Writing to: ' + output_bed_fpath)

    # process configuration
    # for k, v in opts.__dict__.items():
    #     if k.endswith('fpath') and verify_file(v, is_critical=True):
    #         opts.__dict__[k] = verify_file(v, k)
    if cnf.output_grch and cnf.output_hg:
        info('you cannot specify --output-hg and --output-grch simultaneously!')
    # if not which(opts.bedtools):
    #     info('bedtools executable not found, please specify correct path (current is %s)! '
    #         'Did you forget to execute "module load bedtools"?' % opts.bedtools)

    # if opts.debug:
    #     info('Configuration: ')
    #     for k, v in opts.__dict__.items():
    #         info('\t' + k + ': ' + str(v))
    info()

    # opts.ensembl_bed_fpath = verify_file(opts.ensembl_bed_fpath or \
    #     ('/ngs/reference_data/genomes/Hsapiens/' + opts.genome + '/bed/Exons/Exons.with_genes.bed'))

    # opts.refseq_bed_fpath = verify_file(opts.refseq_bed_fpath or \
    #     ('/ngs/reference_data/genomes/Hsapiens/' + opts.genome + '/bed/Exons/RefSeq/RefSeq_CDS_miRNA.all_features.bed'))

    return input_bed_fpath, output_bed_fpath, work_dirpath, cnf


class BedParams:
    GRCh_to_hg = {'MT': 'chrM', 'X': 'chrX', 'Y': 'chrY'}
    for i in range(1, 23):
        GRCh_to_hg[str(i)] = 'chr' + str(i)
    hg_to_GRCh = dict((v, k) for k, v in GRCh_to_hg.items())

    def __init__(self, header=list(), controls=list(), GRCh_names=None, n_cols_needed=None):
        self.header = header
        self.controls = controls
        self.GRCh_names = GRCh_names
        self.n_cols_needed = n_cols_needed

    @staticmethod
    def calc_n_cols_needed(line):
        columns = line.strip().split('\t')
        if len(columns) < 8:
            return 4
        elif columns[6].isdigit() and columns[7].isdigit():  # primers info
            return 8
        else:
            return 4


class Region(SortableByChrom):
    GRCh_names = None
    n_cols_needed = 4
    approved_genes = []
    key_genes = []

    def __init__(self, chrom, chrom_ref_order, start, end):
        SortableByChrom.__init__(self, chrom, chrom_ref_order)
        self.start = start
        self.end = end
        self.symbol = None
        self.type = None
        self.rest = None

    def set_symbol(self, symbol):
        self.symbol = symbol
        if self.symbol in self.key_genes:
            self.type = 'key'
        elif self.symbol in self.approved_genes:
            self.type = 'approved'
        else:
            self.type = 'not_approved'

    def __str__(self):
        fs = [BedParams.hg_to_GRCh[self.chrom] if self.GRCh_names else self.chrom,
              str(self.start), str(self.end), self.symbol] + (self.rest if self.n_cols_needed > 4 else [])
        return '\t'.join(fs)

    def get_key(self):
        return SortableByChrom.get_key(self), self.start, self.end, self.is_control(), self.symbol

    def is_control(self):
        if self.symbol and self.symbol.startswith('_CONTROL'):
            return True
        return False

    def __lt__(self, other):
        return self.get_key() < other.get_key()

    def __eq__(self, other):
        return self.get_key() == other.get_key()

    def __hash__(self):
        return hash(self.get_key())


def _preprocess(cnf, bed_fpath, work_dirpath, chrom_order):
    bed_params = BedParams()
    output_fpath = __intermediate_fname(work_dirpath, bed_fpath, 'prep')
    info('preprocessing: ' + bed_fpath + ' --> ' + output_fpath)
    with open(bed_fpath, 'r') as in_f:
        with open(output_fpath, 'w') as out_f:
            for line in in_f:
                if line.startswith('#') or line.startswith('track') or line.startswith('browser'):  # header
                    bed_params.header.append(line if line.startswith('#') else '#' + line)
                else:
                    cur_ncn = BedParams.calc_n_cols_needed(line)
                    if bed_params.n_cols_needed is not None and cur_ncn != bed_params.n_cols_needed:
                        critical('number and type of columns should be the same on all lines!')
                    bed_params.n_cols_needed = cur_ncn
                    if line.startswith('chr'):
                        if bed_params.GRCh_names is not None and bed_params.GRCh_names:
                            critical('mixing of GRCh and hg chromosome names!')
                        bed_params.GRCh_names = False
                        if line.startswith('chrMT'):  # common misprint, correcting chrMT --> chrM
                            processed_line = '\t'.join(['chrM'] + line.split('\t')[1:])
                        else:
                            processed_line = line
                    elif line.split('\t')[0] in BedParams.GRCh_to_hg:  # GRCh chr names
                        if bed_params.GRCh_names is not None and not bed_params.GRCh_names:
                            critical('mixing of GRCh and hg chromosome names!')
                        bed_params.GRCh_names = True
                        processed_line = '\t'.join([BedParams.GRCh_to_hg[line.split('\t')[0]]] + line.split('\t')[1:])
                    else:
                        critical('incorrect chromosome name!')

                    entries = processed_line.strip().split('\t')
                    chrom = entries[0]
                    start = int(entries[1])
                    end = int(entries[2])
                    r = Region(chrom, chrom_order.get(chrom), start, end)
                    if r.is_control():
                        r.set_symbol(entries[3] if len(entries) > 3 else '{0}:{1}-{2}'.format(chrom, start, end))
                        r.rest = entries[4:] if len(entries) > 4 else None
                        bed_params.controls.append(r)
                    else:
                        out_f.write(processed_line)
    return output_fpath, bed_params


def _annotate(bed_fpath, work_dirpath, cnf):
    annotated_files = []
    input_fpath = bed_fpath
    references = [('RefSeq', cnf.genome.features), ('Ensembl', cnf.genome.ensembl)]

    for id, (db_name, db_bed_fpath) in enumerate(references):
        output_fpath = __intermediate_fname(work_dirpath, bed_fpath, 'ann_' + db_name.lower())
        info('annotating based on {db_name}: {bed_fpath} --> {output_fpath}'.format(**locals()))
        annotate_bed_py = sys.executable + ' ' + splitext(annotate_bed.__file__)[0] + '.py'

        cmdline = '{annotate_bed_py} {input_fpath} --reference {db_bed_fpath} -o {output_fpath} --genome {cnf.genome}'.format(**locals())
        __call(cnf, cmdline)

        if id < len(references) - 1:
            if cnf.debug:
                info("filtering annotated and not annotated regions into separate files:")
            only_annotated_bed = __intermediate_fname(work_dirpath, bed_fpath, 'only_ann_' + db_name.lower())
            not_annotated_bed = __intermediate_fname(work_dirpath, bed_fpath, 'not_ann_' + db_name.lower())
            with open(only_annotated_bed, 'w') as out:
                cmdline = 'grep -v -E "\.$" {output_fpath}'.format(**locals())
                if cnf.debug:
                    info(cmdline + ' > ' + only_annotated_bed)
                subprocess.call(cmdline, shell=True, stdout=out)
            with open(not_annotated_bed, 'w') as out:
                cmdline = 'grep -E "\.$" {output_fpath}'.format(**locals())
                if cnf.debug:
                    info(cmdline + ' > ' + not_annotated_bed)
                subprocess.call(cmdline, shell=True, stdout=out)
            if not cnf.debug:
                os.remove(output_fpath)
            output_fpath = only_annotated_bed
            input_fpath = not_annotated_bed
        annotated_files.append(output_fpath)
        if id != 0 and not cnf.debug:
            os.remove(input_fpath)

    return annotated_files


def _postprocess(input_fpath, annotated_fpaths, bed_params, output_bed_fpath, cnf, chrom_order):
    '''
    1. Sorts.
    1. Chooses appropriate number of columns (4 or 8 for BEDs with primers).
    2. Removes duplicates.
    '''
    info('postprocessing (sorting, cutting, removing duplicates)')

    key_genes = []
    with open(adjust_path(cnf.key_genes), 'r') as f:
        for line in f:
            key_genes.append(line.strip())
    approved_genes = []
    if cnf.hgnc:
        with open(adjust_path(cnf.hgnc), 'r') as f:
            f.readline()  # header
            for line in f:
                approved_genes.append(line.split('\t')[0])

    Region.GRCh_names = bed_params.GRCh_names
    if cnf.output_grch:
        Region.GRCh_names = True
        if cnf.debug and not bed_params.GRCh_names:
            info('Changing chromosome names from hg-style to GRCh-style.')
    if cnf.output_hg:
        Region.GRCh_names = False
        if cnf.debug and bed_params.GRCh_names:
            info('Changing chromosome names from GRCh-style to hg-style.')
    Region.n_cols_needed = bed_params.n_cols_needed
    Region.key_genes = key_genes
    Region.approved_genes = approved_genes

    input_regions = set()  # we want only unique regions
    with open(adjust_path(input_fpath)) as f:
        for line in f:
            entries = line.strip().split('\t')
            chrom = entries[0]
            start = int(entries[1])
            end = int(entries[2])
            r = Region(chrom, chrom_order.get(chrom), start, end)
            r.set_symbol(entries[3] if len(entries) > 3 else '{0}:{1}-{2}'.format(chrom, start, end))
            r.rest = entries[4:] if len(entries) > 4 else None
            input_regions.add(r)

    annotated_regions = []
    for annotated_fpath in annotated_fpaths:
        with open(adjust_path(annotated_fpath)) as f:
            for line in f:
                entries = line.strip().split('\t')
                chrom = entries[0]
                start = int(entries[1])
                end = int(entries[2])
                r = Region(chrom, chrom_order.get(chrom), start, end)
                r.set_symbol(entries[3] if len(entries) > 3 else '{0}:{1}-{2}'.format(chrom, start, end))
                r.rest = entries[4:] if len(entries) > 4 else None
                annotated_regions.append(r)

    # starting to output result
    with open(adjust_path(output_bed_fpath), 'w') as f:
        for line in bed_params.header:
            f.write(line)

        annotated_regions.sort()
        i = 0
        prev_region = None
        not_a_gene_count = 0
        solid_regions = []
        prev_is_solid = False
        all_regions = []
        for cur_region in sorted(list(input_regions) + bed_params.controls):
            if not cur_region.is_control():
                assert annotated_regions[i] == cur_region, str(cur_region) + ' != ' + str(annotated_regions[i]) + '(i=%d)' % i
                if annotated_regions[i].symbol != '.':
                    cur_region.set_symbol(annotated_regions[i].symbol)
                else:
                    if prev_region is None or \
                       prev_region.chrom != cur_region.chrom or not prev_region.symbol.startswith("not_a_gene"):
                        not_a_gene_count += 1
                    cur_region.set_symbol("not_a_gene_%d" % not_a_gene_count)
                i += 1
                ambiguous_regions = [cur_region]
                while i < len(annotated_regions) and annotated_regions[i] == cur_region:  # processing duplicates
                    if annotated_regions[i].symbol != '.' and annotated_regions[i].symbol != cur_region.symbol:
                        duplicate = copy.deepcopy(cur_region)
                        duplicate.set_symbol(annotated_regions[i].symbol)
                        if duplicate.type == 'approved' and cur_region.type == 'not_approved':
                            cur_region = duplicate
                            ambiguous_regions = [cur_region]
                        elif annotated_regions[i].type == 'key' and cur_region.type != 'key':
                            cur_region = duplicate
                            ambiguous_regions = [cur_region]
                            if cnf.debug:
                                info('key gene priority over approved gene was used')
                        elif annotated_regions[i].type == cur_region.type:
                            ambiguous_regions.append(duplicate)
                    i += 1
                if len(ambiguous_regions) == 1:
                    if not prev_is_solid:
                        solid_regions.append(cur_region)
                    prev_is_solid = True
                    all_regions.append(cur_region)
                else:
                    if prev_is_solid:
                        solid_regions.append(prev_region)
                    prev_is_solid = False
                    all_regions.append(ambiguous_regions)
            else:
                all_regions.append(cur_region)
            prev_region = cur_region

        # outputting results
        cur_solid_id = -1
        for entry in all_regions:
            if isinstance(entry, list):  # list of ambiguous regions
                cur_region = entry[0]
                while cur_solid_id + 1 < len(solid_regions) and cur_region > solid_regions[cur_solid_id + 1]:
                    cur_solid_id += 1
                found = False
                if cur_solid_id >= 0 and cur_region > solid_regions[cur_solid_id] \
                        and cur_region.chrom == solid_regions[cur_solid_id].chrom:
                    prev_solid = solid_regions[cur_solid_id]
                    for cur_region in entry:
                        if cur_region.symbol == prev_solid.symbol:
                            found = True
                            if cnf.debug:
                                info('gene name was chosen based on previous solid region')
                            break
                if not found and cur_solid_id + 1 < len(solid_regions) and cur_region < solid_regions[cur_solid_id + 1] \
                        and cur_region.chrom == solid_regions[cur_solid_id + 1].chrom:
                    next_solid = solid_regions[cur_solid_id + 1]
                    for cur_region in entry:
                        if cur_region.symbol == next_solid.symbol:
                            found = True
                            if cnf.debug:
                                info('gene name was chosen based on next solid region')
                            break
                if not found:
                    cur_region = entry[0]
            else:
                cur_region = entry
            f.write(str(cur_region) + '\n')  # automatically outputs correct number of columns and GRCh/hg names


def __intermediate_fname(work_dir, fname, suf):
    output_fname = add_suffix(fname, suf)
    return join(work_dir, basename(output_fname))


def __call(cnf, cmdline, output_fpath=None):
    stdout = open(output_fpath, 'w') if output_fpath else None
    stderr = None if cnf.debug else open('/dev/null', 'w')
    if cnf.debug:
        info(cmdline)
    ret_code = subprocess.call(cmdline, shell=True, stdout=stdout, stderr=stderr, stdin=None)
    return ret_code


def main():
    input_bed_fpath, output_bed_fpath, work_dirpath, cnf = _read_args(sys.argv[1:])

    chr_lengths = get_chr_lengths_from_seq(cnf.genome.seq)
    chrom_order = {c: i for i, (c, l) in enumerate(chr_lengths)}

    preprocessed_fpath, bed_params = _preprocess(cnf, input_bed_fpath, work_dirpath, chrom_order)
    annotated_fpaths = _annotate(preprocessed_fpath, work_dirpath, cnf)

    _postprocess(preprocessed_fpath, annotated_fpaths, bed_params, output_bed_fpath, cnf, chrom_order)
    if not cnf.debug:
        for f in [preprocessed_fpath] + annotated_fpaths:
            os.remove(f)

    try:
        shutil.rmtree(work_dirpath)
    except OSError:
        pass

if __name__ == '__main__':
    main()
