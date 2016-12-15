import os
import sys
import math
from collections import OrderedDict
from os.path import isfile, join, abspath, basename, dirname, getctime, getmtime, splitext, realpath
from subprocess import check_output

from ngs_utils.bedtools import BedTool
from ngs_utils import call_process
from ngs_utils.file_utils import intermediate_fname, iterate_file, splitext_plus, verify_file, adjust_path, add_suffix, \
    safe_mkdir, file_transaction, which, file_exists, open_gzipsafe, can_reuse
from ngs_utils.logger import info, critical, warn, err, debug
from ngs_utils import reference_data as ref
from ngs_utils.utils import md5


def get_chrom_order(genome=None, fai_fpath=None):
    chr_lengths = ref.get_chrom_lengths(genome, fai_fpath)
    chr_order = {c: i for i, (c, l) in enumerate(chr_lengths)}
    return chr_order


class SortableByChrom:
    def __init__(self, chrom, chrom_ref_order):
        self.chrom = chrom
        if isinstance(chrom_ref_order, dict):
            self.chrom_ref_order = chrom_ref_order.get(chrom, -1)
        else:
            self.chrom_ref_order = chrom_ref_order

    def get_key(self):
        return self.chrom_ref_order


class Region(SortableByChrom):
    def __init__(self, chrom, start, end, other_fields, chrom_ref_order):
        SortableByChrom.__init__(self, chrom, chrom_ref_order)
        self.start = int(start)
        self.end = int(end) if end is not None else None
        self.other_fields = list(other_fields)

    def get_key(self):
        return self.chrom_ref_order, self.start, self.end, self.other_fields

    def __str__(self):
        ts = [self.chrom, self.start, self.end] + self.other_fields
        return '"' + '\t'.join(map(str, ts)) + '"'

    def __repr__(self):
        return self.__str__()


def bam_to_bed(cnf, bam_fpath, to_gzip=True):
    info('Converting the BAM to BED to save some memory.')  # from here: http://davetang.org/muse/2015/08/05/creating-a-coverage-plot-using-bedtools-and-r/
    bam_bed_fpath = splitext_plus(bam_fpath)[0] + ('.bed.gz' if to_gzip else '.bed')
    bedtools = which('bedtools')
    gzip = which('gzip')
    cmdline = '{bedtools} bamtobed -i {bam_fpath}'.format(**locals())
    cmdline += ' | {gzip}'.format(**locals()) if to_gzip else ''
    call_process.run(cmdline, output_fpath=bam_bed_fpath)
    return bam_bed_fpath


def bam_to_bed_nocnf(bam_fpath, bedtools='bedtools', gzip='gzip'):
    info('Converting the BAM to BED to save some memory.')  # from here: http://davetang.org/muse/2015/08/05/creating-a-coverage-plot-using-bedtools-and-r/
    bam_bed_fpath = splitext_plus(bam_fpath)[0] + '.bed.gz'
    cmdline = '{bedtools} bamtobed -i {bam_fpath} | {gzip} > {bam_bed_fpath}'.format(**locals())
    info(cmdline)
    os.system(cmdline)
    bam_bed_fpath = verify_file(bam_bed_fpath)
    if bam_bed_fpath:
        info('Done, saved to ' + bam_bed_fpath)
    else:
        err('Error, result is non-existent or empty')
    return bam_bed_fpath


def count_bed_cols(bed_fpath):
    with open(bed_fpath) as f:
        for l in f:
            if l and l.strip() and not l.startswith('#'):
                return len(l.split('\t'))
    # return len(next(dropwhile(lambda x: x.strip().startswith('#'), open(bed_fpath))).split('\t'))
    err('Empty bed file: ' + bed_fpath)
    return None


def merge_overlaps(work_dir, bed_fpath, distance=None):
    """Merge bed file intervals to avoid overlapping regions.
    Overlapping regions (1:1-100, 1:90-100) cause issues with callers like FreeBayes
    that don't collapse BEDs prior to using them.
    """
    output_fpath = intermediate_fname(work_dir, bed_fpath, 'merged')
    if isfile(output_fpath) and verify_file(output_fpath, cmp_f=bed_fpath):
        return output_fpath

    with file_transaction(work_dir, output_fpath) as tx:
        kwargs = dict(d=distance) if distance else dict()
        BedTool(bed_fpath).merge(**kwargs).saveas(tx)
    return output_fpath


def remove_comments(work_dir, bed_fpath, reuse=False):
    def f(l, i):
        if not l.startswith('#'):
            return l
        else:
            return None
    return iterate_file(work_dir, bed_fpath, f, suffix='rmcmt', reuse=reuse)


def calc_region_number(bed_fpath):
    with open(bed_fpath) as f:
        return sum(1 for l in f if l.strip() and not l.strip().startswith('#'))


def get_genes_from_bed(bed_fpath, chrom_index=0, gene_index=3):
    gene_keys_set = set()
    gene_keys_list = list()
    with open(bed_fpath) as f:
        for line in f:
            if not line or not line.strip() or line.startswith('#'):
                continue

            tokens = line.replace('\n', '').split('\t')
            if len(tokens) <= gene_index or len(tokens) <= chrom_index:
                continue

            chrom = tokens[chrom_index]
            for gn in tokens[gene_index].split(','):
                if (gn, chrom) not in gene_keys_set:
                    gene_keys_set.add((gn, chrom))
                    gene_keys_list.append((gn, chrom))

    return gene_keys_set, gene_keys_list


# def group_and_merge_regions_by_gene(work_dir, bed_fpath, keep_genes=False):
#     output_fpath = intermediate_fname(work_dir, bed_fpath, 'grp_mrg')
#
#     group_merge_bed_py = join(tc.code_base_path, 'tools', 'bed_processing', 'group_and_merge_by_gene.py')
#
#     cmdline = '{group_merge_bed_py} {bed_fpath}'.format(**locals())
#     if not keep_genes:
#         cmdline += ' | grep -vw Gene'
#
#     call_process.run(cmdline, 'Group and merge bed', output_fpath)
#
#     return output_fpath


def cut(fpath, col_num, output_fpath=None, reuse=False):
    output_fpath = output_fpath or add_suffix(fpath, 'cut')
    cmdline = 'cut -f' + ','.join(map(str, range(1, col_num + 1))) + ' ' + fpath + ' > ' + output_fpath
    call_process.run(cmdline, output_fpath=output_fpath, reuse=reuse)
    return output_fpath


def bgzip_and_tabix(fpath, reuse=False, tabix_parameters='', **kwargs):
    gzipped_fpath = join(fpath + '.gz')
    tbi_fpath = gzipped_fpath + '.tbi'

    if reuse and \
           file_exists(gzipped_fpath) and (getctime(gzipped_fpath) >= getctime(fpath) if file_exists(fpath) else True) and \
           file_exists(tbi_fpath) and getctime(tbi_fpath) >= getctime(gzipped_fpath):
        info('Actual compressed file and index exist, reusing')
        return gzipped_fpath

    info('Compressing and tabixing file, writing ' + gzipped_fpath + '(.tbi)')
    bgzip = which('bgzip')
    tabix = which('tabix')
    if not bgzip:
        err('Cannot index file because bgzip is not found')
    if not tabix:
        err('Cannot index file because tabix is not found')
    if not bgzip and not tabix:
        return fpath

    if isfile(gzipped_fpath):
        os.remove(gzipped_fpath)
    if isfile(tbi_fpath):
        os.remove(tbi_fpath)

    info('BGzipping ' + fpath)
    cmdline = '{bgzip} {fpath}'.format(**locals())
    call_process.run(cmdline)

    info('Tabixing ' + gzipped_fpath)
    cmdline = '{tabix} {tabix_parameters} {gzipped_fpath}'.format(**locals())
    call_process.run(cmdline)

    return gzipped_fpath


def prep_bed_for_seq2c(work_dir, bed_fpath, reuse=False):
    info()
    info('Doing some Seq2C specific preparation of the BED file...')

    cols = count_bed_cols(bed_fpath)

    seq2c_bed = None
    if 8 > cols > 4:
        seq2c_bed = cut(bed_fpath, 4, reuse=False)
    elif cols > 8:
        seq2c_bed = cut(bed_fpath, 8, reuse=False)
    else:
        seq2c_bed = bed_fpath

    if cols >= 4:
        # removing regions with no gene annotation
        def f(l, i):
            if l.split('\t')[3].strip() == '.': return None
            else: return l
        seq2c_bed = iterate_file(work_dir, seq2c_bed, f, suffix='filt', reuse=reuse)

    info('Done: ' + seq2c_bed)
    return seq2c_bed


def filter_bed_with_gene_set(work_dir, bed_fpath, gene_keys_set, suffix=None, reuse=False):
    met_genes = set()

    def fun(r):
        if not r.name:
            return None
        new_gns = []
        for g in r.name.split(','):
            if (g, r.chrom) in gene_keys_set:
                new_gns.append(g)
                met_genes.add((g, r.chrom))
        if new_gns:
            r.name = ','.join(new_gns)

    return BedTool(bed_fpath).each(fun).saveas().fn, met_genes

    # res_fpath = iterate_file(work_dir, bed, fn, suffix=suffix or 'filt_genes', check_result=False)


def sort_bed(input_bed_fpath, output_bed_fpath=None, work_dir=None, fai_fpath=None, chr_order=None, genome=None, reuse=False):
    input_bed_fpath = verify_bed(input_bed_fpath, is_critical=True)
    output_bed_fpath = adjust_path(output_bed_fpath) if output_bed_fpath \
        else intermediate_fname(work_dir, input_bed_fpath, 'sorted')

    regions = []

    if not chr_order:
        if fai_fpath:
            fai_fpath = verify_file(fai_fpath)
        elif genome:
            fai_fpath = verify_file(ref.get_fai(genome))
        else:
            critical('Either of chr_order, fai_fpath, or genome build name must be specified')
        chr_order = get_chrom_order(fai_fpath=fai_fpath)

    debug('Sorting regions in ' + str(input_bed_fpath))
    if reuse and isfile(output_bed_fpath) and verify_file(output_bed_fpath, cmp_f=input_bed_fpath):
        debug(output_bed_fpath + ' exists, reusing')
        return output_bed_fpath

    with open(input_bed_fpath) as f:
        with file_transaction(work_dir, output_bed_fpath) as tx:
            with open(tx, 'w') as out:
                for l in f:
                    if not l.strip():
                        continue
                    if l.strip().startswith('#'):
                        out.write(l)
                        continue

                    fs = l.strip().split('\t')
                    chrom = fs[0]
                    start = int(fs[1])
                    end = int(fs[2])
                    other_fields = fs[3:]
                    order = chr_order.get(chrom, -1)
                    regions.append(Region(chrom, start, end, other_fields, order))

                for region in sorted(regions, key=lambda r: r.get_key()):
                    fs = [region.chrom, str(region.start), str(region.end)]
                    fs.extend(region.other_fields)
                    out.write('\t'.join(fs) + '\n')

    debug('Sorted ' + str(len(regions)) + ' regions, saved to ' + output_bed_fpath)
    return output_bed_fpath


def calc_sum_of_regions(bed_fpath):
    total_bed_size = 0

    with open(bed_fpath) as f:
        for l in f:
            l = l.strip()
            if not l.startswith('#'):
                start, end = [int(f) for f in l.split('\t')[1:3]]
                total_bed_size += end - start

    return total_bed_size


def get_total_bed_size(bed_fpath):
    return sum(len(x) for x in BedTool(bed_fpath).merge())


def bedtools_version(bedtools):
    v = check_output([bedtools, '--version'])  # bedtools v2.24.0
    try:
        v = int(v.split(' ')[1].split('.')[1])
    except:
        return None
    else:
        return v


this_py_fpath = splitext(__file__)[0] + '.py'
this_py_real_fpath = realpath(this_py_fpath)
project_dir = dirname(this_py_real_fpath)
code_base_path = abspath(join(abspath(project_dir)))



def get_padded_bed_file(work_dir, bed, padding, fai_fpath):
    genome_fpath = fai_fpath
    info('Making bed file for padded regions...')
    bedtools = which('bedtools')
    cmdline = '{bedtools} slop -i {bed} -g {genome_fpath} -b {padding}'.format(**locals())
    output_fpath = intermediate_fname(work_dir, bed, 'padded')
    call_process.run(cmdline, output_fpath=output_fpath)
    return output_fpath


def intersect_bed(work_dir, bed1, bed2):
    bed1_fname, _ = splitext_plus(basename(bed1))
    bed2_fname, _ = splitext_plus(basename(bed2))
    output_fpath = join(work_dir, bed1_fname + '__' + bed2_fname + '.bed')
    bedtools = which('bedtools')
    cmdline = '{bedtools} intersect -u -a {bed1} -b {bed2}'.format(**locals())
    call_process.run(cmdline, output_fpath=output_fpath, checks=[call_process.file_exists])
    return output_fpath


def verify_bed(bed, description='', is_critical=False, silent=False):
    if isinstance(bed, BedTool):
        return bed

    fpath = adjust_path(bed)
    if not verify_file(fpath, description, is_critical=is_critical, silent=silent):
        return None

    error = BedFile(fpath).checkformat()
    if error:
        fn = critical if is_critical else err
        fn('Error: incorrect bed file format (' + fpath + '): ' + str(error) + '\n')
        return None

    return fpath


def clean_bed(bed_fpath, work_dir):
    clean_fpath = intermediate_fname(work_dir, bed_fpath, 'clean')

    if not can_reuse(clean_fpath, bed_fpath):
        bed = BedTool(bed_fpath)
        bed = bed.filter(lambda x: x.chrom and
                         not any(x.chrom.startswith(e) for e in ['#', ' ', 'track', 'browser']))
        bed = bed.remove_invalid()
        with file_transaction(work_dir, clean_fpath) as tx_out_file:
            bed.saveas(tx_out_file)
        verify_bed(clean_fpath, is_critical=True)
        debug('Saved Seq2C clean BED file into ' + clean_fpath)
    return clean_fpath


def check_md5(work_dir, fpath, file_type, silent=False):
    md5_fpath = join(work_dir, file_type + '_md5.txt')
    new_md5 = md5(fpath)
    info('md5 of ' + fpath + ' is ' + str(new_md5))
    prev_md5 = None
    if isfile(md5_fpath):
        with open(md5_fpath) as f:
            prev_md5 = f.read()
    else:
        info('Previous md5 file ' + md5_fpath + ' does not exist')
    info('Checking previous md5 from ' + md5_fpath + ': ' + str(prev_md5))

    if prev_md5 == new_md5:
        if not silent:
            debug('Reusing previous ' + file_type.upper() + ' files.')
        return True
    else:
        if not silent:
            info('Pre-processing input ' + file_type.upper() + ' file')
        if prev_md5:
            if not silent:
                info('Prev ' + file_type.upper() + ' md5: ' + str(prev_md5))
                info('New ' + file_type.upper() + ' md5: ' + str(new_md5))

        with open(md5_fpath, 'w') as f:
            f.write(str(new_md5))
        return False


class BedFile:
    def __init__(self, _filename):
        self.filename = _filename
        self.chrs = None
        self.nregions = None

    def checkformat(self):
        """************************************************************************************************************************************************************
        Task: checks the format of the bed file. The only requirements checked are that each line presents at least 3 tab separated columns, the
            two on the right must present integer values indicating the start/end position respectively. Right value must be greater than the
            left value.
        Outputs:
            err: string containing the detected error. Empty string in case of a correct format.
        ************************************************************************************************************************************************************"""

        fd = open_gzipsafe(self.filename)

        line = fd.readline()
        while line.startswith('#'):
            line = fd.readline()

        fields = line.split('\t')
        lc = 1
        error = ''

        # Checks that the two columns on the right contain integer values
        try:
            # Parses each line and checks that there are at least 3 fields, the two on the right containing integer values and being the right one
            # greater than the left one
            while line <> '' and len(fields) > 2 and int(fields[1]) <= int(fields[2]):
                lc += 1
                line = fd.readline()
                fields = line.split('\t')
        except ValueError:
            error += 'Incorrect start/end values at line ' + str(lc) + '\n'
            error += 'Start/End coordinates must be indicated with integer values. The right value must be greater than the left value.\n'
            error += 'Line found: ' + line
            fd.close()

            return error

        # If it get to this point means that either the file ended or there is a line with less than 3 fields
        if line <> '':
            error += 'Incorrect line format at line ' + str(lc) + '\n'
            error += 'At least three columns are expected in each line\n'
            error += 'The right value must be greater than the left value.\n'
            error += 'Line found: ' + line
            fd.close()

        return error

    def count_lines(self, filename=None):
        if filename is None:
            filename = self.filename
        return len(open(filename).readlines())


# class Region:
#     def __init__(self, sample_name=None, gene_name=None, transcript_id=None,
#                  exon_num=None, strand=None, biotype=None, feature=None, extra_fields=list(),
#                  chrom=None, start=None, end=None, size=None, min_depth=None, median_depth=None,
#                  avg_depth=None, std_dev=None, rate_within_normal=None, bases_by_depth=None):
#
#         self.sample_name = sample_name
#         self.gene_name = gene_name
#         self.feature = feature
#         self.extra_fields = extra_fields  # for exons, extra_fields is [Gene, Exon number, Strand]
#         self.exon_num = exon_num
#         self.strand = strand
#         self.biotype = biotype
#         self.transcript_id = transcript_id
#
#         self.chrom = chrom
#         self.start = start  # int
#         self.end = end      # int
#         self.size = size    # int
#         self.min_depth = min_depth
#         self.bases_by_depth = bases_by_depth or defaultdict(int)  # filled in from the "bedcoverage hist" output
#
#         # Calculated once on "sum_up()", when all self.bases_by_depth are there:
#         self.avg_depth = avg_depth  # float
#         self.median_depth = median_depth  # float
#         self.std_dev = std_dev
#         self.rate_within_normal = rate_within_normal
#         self.bases_within_threshs = None    # OrderedDict((depth, 0) for depth in depth_thresholds)
#         self.rates_within_threshs = None    # defaultdict(float)
#
#         self.missed_by_db = dict()
#         self.var_num = None
#
#     # def get_chrom_num(self):
#     #     digits = [c for c in self.chrom if c.isdigit()]
#     #     if digits:
#     #         return int(''.join(digits))
#     #     if 'M' in self.chrom:
#     #         return 0
#     #     if 'X' in self.chrom:
#     #         return 23
#     #     if 'Y' in self.chrom:
#     #         return 24
#     #     else:
#     #         return 25
#     #
#     # @staticmethod
#     # def get_order_key(r):
#     #     return r.get_chrom_num(), r.get_start(), r.get_end(), r.gene_name
#
#     def get_start(self):
#         return self.start
#
#     def get_end(self):
#         return self.end
#
#     def get_avg_depth(self):
#         return self.avg_depth
#
#     def get_std_dev(self):
#         return self.std_dev
#
#     def get_bases_by_depth(self):
#         return self.bases_by_depth
#
#     def get_size(self):
#         if self.size is not None:
#             return self.size
#         if self.start and self.end:
#             return abs(self.end - self.start)
#         return None
#
#     def add_bases_for_depth(self, depth, bases):
#         if depth in self.bases_by_depth:
#             err('Duplicated depth ' + str(depth) + ' for the region ' + str(self))
#         else:
#             self.bases_by_depth[depth] += bases
#
#     def __str__(self):
#         ts = [self.sample_name, self.chrom, self.start, self.end,
#               self.gene_name, self.feature]
#         return '"' + '\t'.join(map(str, ts)) + '"'
#
#     def to_bed_str(self):
#         return '\t'.join([self.chrom, str(self.start), str(self.end), self.gene_name] + '\n')
#
#     def __repr__(self):
#         return self.__str__()
#
#     def calc_bases_within_threshs(self, depth_thresholds):
#         if self.bases_within_threshs is not None:
#             return self.bases_within_threshs
#
#         if self.bases_by_depth is None:
#             err('Error: self.bases_by_depth is None for ' + str(self))
#
#         self.bases_within_threshs, self.rates_within_threshs = calc_bases_within_threshs(
#             self.bases_by_depth, self.get_size(), depth_thresholds)
#
#         return self.bases_within_threshs
#
#     def calc_avg_depth(self):
#         if self.avg_depth is not None:
#             return self.avg_depth
#
#         if self.bases_by_depth:
#             depth_sum = sum(
#                 depth * bases
#                 for depth, bases
#                 in self.bases_by_depth.items())
#             self.avg_depth = float(depth_sum) / self.get_size() if self.get_size() else None
#             return self.avg_depth
#
#     def calc_std_dev(self, avg_depth):
#         if avg_depth is None:
#             return None
#
#         if self.std_dev is not None:
#             return self.std_dev
#
#         if self.bases_by_depth:
#             sum_of_sq_var = sum(
#                 (depth - avg_depth) ** 2 * bases
#                 for depth, bases
#                 in self.bases_by_depth.items())
#             sz = self.get_size()
#             if sz and sz > 0:
#                 d = float(sum_of_sq_var) / float(sz)
#                 try:
#                     self.std_dev = math.sqrt(d)
#                 except ValueError, e:
#                     print 'float(sum_of_sq_var) =', float(sum_of_sq_var)
#                     print 'float(sz) =', float(sz)
#                     print 'd =', d
#                     print self.sample_name, self.gene_name, self.chrom, ':', self.start, '-', self.end
#                     # print 'math.sqrt(d) =', math.sqrt(d)
#                     critical(str(e))
#             return self.std_dev
#
#     def calc_rate_within_normal(self, avg_depth):
#         if avg_depth is None:
#             return None
#
#         if self.rate_within_normal is not None:
#             return self.rate_within_normal
#
#         return calc_rate_within_normal(self.bases_by_depth, avg_depth, self.get_size())
#
#     def sum_up(self, depth_thresholds):
#         self.calc_avg_depth()
#         self.calc_std_dev(self.avg_depth)
#         self.calc_bases_within_threshs(depth_thresholds)
#         self.calc_rate_within_normal(self.avg_depth)
#         return self.bases_within_threshs,\
#                self.bases_by_depth, \
#                self.avg_depth, \
#                self.std_dev, \
#                self.rate_within_normal
#
#     def intersect(self, reg2):
#         return self.chrom == reg2.chrom and \
#                (reg2.start < self.start < reg2.end or
#                 self.start < reg2.start < self.end)


def calc_rate_within_normal(bases_by_depth, avg_depth, total_size):
    bases_within_normal = sum(
        bases
        for depth, bases
        in bases_by_depth.items()
        if math.fabs(avg_depth - depth) < 0.2 * avg_depth)

    return 1.0 * bases_within_normal / total_size if total_size else None


def calc_bases_within_threshs(bases_by_depth, total_size, depth_thresholds):
    bases_within_threshs = OrderedDict((depth, 0) for depth in depth_thresholds)
    rates_within_threshs = OrderedDict((depth, None) for depth in depth_thresholds)

    for depth, bases in bases_by_depth.iteritems():
        for t in depth_thresholds:
            if depth >= t:
                bases_within_threshs[t] += bases
    for t in depth_thresholds:
        bs = bases_within_threshs[t]
        if total_size > 0:
            rate = 1.0 * bases_within_threshs[t] / total_size
            if rate > 1:
                critical('Error: rate is > 1: rate = ' + str(rate) + ', bases = ' + str(bs) + ', size = ' + str(total_size))
            rates_within_threshs[t] = rate

    return bases_within_threshs, rates_within_threshs


class GeneInfo(Region):
    """ Collects assisiated exons and overlapping amlicons.
        - Knows its sample, gene name and chromosome.

        - Stores amplicons and exons ("Region" instances).

        - Supports extending with exons in sorted by starting position order;
          when adding a new exon, recalculates start, end, size and based_by_depth.
    """
    def __init__(self, gene_name, chrom=None, strand=None, feature='Gene-Exon', exon_num=None):
        Region.__init__(self,
            sample_name=None, gene_name=gene_name, exon_num=exon_num, strand=strand,
            feature=feature, chrom=chrom)
        self.exons = []
        self.amplicons = []
        self.non_overlapping_exons = []
        self.size = 0
        self.min_depth = None
        # self.amplicon_gene_info = GeneInfo(
        #     sample_name=sample_name, gene_name=gene_name, exon_num=exon_num, strand=strand,
        #     feature='Gene-Capture', chrom=chrom)

    def get_exons(self):
        return self.exons  # self.subregions_by_feature['Exon']['regions']

    def get_amplicons(self):
        return self.amplicons  # self.subregions_by_feature['Capture']['regions']

    def add_exon(self, exon):  # exons come sorted by start
        # Start, End and Size are set earlier when parsing the features file
        # if self.end and self.end > exon.start:
        #     self.size += exon.end - self.end
        # else:
        #     self.size += exon.get_size()
        # self.start = min(self.start, exon.start) if self.start else exon.start
        # self.end = max(self.end, exon.end) if self.end else exon.end
        self.exons.append(exon)
        for depth, bases in exon.bases_by_depth.items():
            self.bases_by_depth[depth] += bases

        self.min_depth = min(self.min_depth, exon.min_depth) if self.min_depth else exon.min_depth

    def add_amplicon(self, amplicon):
        # amplicon = copy.copy(amplicon)
        amplicon.gene_name = amplicon.gene_name or self.gene_name
        self.amplicons.append(amplicon)
        if self.gene_name == '.':  # no exons!
            if self.end and self.end > amplicon.start:
                self.size += amplicon.end - self.end
            else:
                self.size += amplicon.get_size()
            self.start = min(self.start, amplicon.start) if self.start else amplicon.start
            self.end = max(self.end, amplicon.end) if self.end else amplicon.end
            for depth, bases in amplicon.bases_by_depth.items():
                self.bases_by_depth[depth] += bases
            self.min_depth = min(self.min_depth, amplicon.min_depth) if self.min_depth else amplicon.min_depth


def build_gene_objects_list(features_bed, gene_keys_list):
    # info('Making unique gene list without affecting the order')
    # fixed_gene_names_list = []
    # added_gene_names_set = set()
    # for i in range(len(gene_names_list)):
    #     gene_name = gene_names_list[i]
    #     if gene_name not in added_gene_names_set:
    #         fixed_gene_names_list.append(gene_name)
    #         added_gene_names_set.add(gene_name)
    # gene_names_list = fixed_gene_names_list
    # info('Uniq gene list contains ' + str(len(gene_names_list)) + ' genes')
    gene_by_name_and_chrom = OrderedDict()

    if gene_keys_list:
        debug()
        debug('Initiate the Gene object dict')
        for gn, chrom in gene_keys_list:
            gene_by_name_and_chrom[(gn, chrom)] = GeneInfo(gene_name=gn, chrom=chrom)
        debug('Processed ' + str(len(gene_keys_list)) + ' gene records -> ' + str(len(gene_by_name_and_chrom)) + ' uniq gene sybmols')

    if features_bed and gene_by_name_and_chrom:
        # info()
        # info('Filtering exon bed file to have only gene records...')
        # exons_only_genes_bed = intermediate_fname(cnf, exons_bed, 'only_genes')
        # call(cnf, 'grep -w Gene ' + exons_bed, output_fpath=exons_only_genes_bed)
        # info('Saved genes to ' + exons_only_genes_bed)

        debug()
        debug('Setting start and end for the genes (based only on the target gene names found in the features list)')
        def fun(c):
            chrom, start, end, symbol, _, strand, feature, biotype, transcript = c.fields
            if feature == 'Transcript':
                gene_by_name_and_chrom[(symbol, chrom)].chrom = chrom
                gene_by_name_and_chrom[(symbol, chrom)].start = int(start)
                gene_by_name_and_chrom[(symbol, chrom)].end = int(end)
                gene_by_name_and_chrom[(symbol, chrom)].size = int(end) - int(start)
                gene_by_name_and_chrom[(symbol, chrom)].strand = strand
                gene_by_name_and_chrom[(symbol, chrom)].biotype = biotype
                gene_by_name_and_chrom[(symbol, chrom)].transcript_id = transcript
        features_bed.each(fun)
        debug()

    if not gene_by_name_and_chrom:
        critical('Err: no genes in the list')

    return gene_by_name_and_chrom


def proc_regions(regions, fn, *args, **kwargs):
    i = 0
    for region in regions:
        i += 1

        fn(region, *args, **kwargs)

        if i % 10000 == 0:
            info('Processed {0:,} regions.'.format(i))

    if not i % 10000:
        info('Processed {0:,} regions.'.format(i))

