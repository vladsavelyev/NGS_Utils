import math
import os
from collections import defaultdict, OrderedDict
from os.path import isfile

from Utils.file_utils import file_transaction, verify_file, adjust_path
from Utils.logger import info, err, critical
import Utils.reference_data as ref


def get_chrom_order(genome=None, fai_fpath=None):
    chr_lengths = ref.get_chrom_lengths(genome, fai_fpath)
    chr_order = {c: i for i, (c, l) in enumerate(chr_lengths)}
    return chr_order


class SortableByChrom:
    def __init__(self, chrom, chrom_ref_order):
        self.chrom = chrom
        self.chrom_ref_order = chrom_ref_order

    def get_key(self):
        return self.chrom_ref_order


class Region:
    def __init__(self, sample_name=None, gene_name=None, transcript_id=None,
                 exon_num=None, strand=None, biotype=None, feature=None, extra_fields=list(),
                 chrom=None, start=None, end=None, size=None, min_depth=None,
                 avg_depth=None, std_dev=None, rate_within_normal=None, bases_by_depth=None):

        self.sample_name = sample_name
        self.gene_name = gene_name
        self.feature = feature
        self.extra_fields = extra_fields  # for exons, extra_fields is [Gene, Exon number, Strand]
        self.exon_num = exon_num
        self.strand = strand
        self.biotype = biotype
        self.transcript_id = transcript_id

        self.chrom = chrom
        self.start = start  # int
        self.end = end      # int
        self.size = size    # int
        self.min_depth = min_depth
        self.bases_by_depth = bases_by_depth or defaultdict(int)  # filled in from the "bedcoverage hist" output

        # Calculated once on "sum_up()", when all self.bases_by_depth are there:
        self.avg_depth = avg_depth  # float
        self.std_dev = std_dev
        self.rate_within_normal = rate_within_normal
        self.bases_within_threshs = None    # OrderedDict((depth, 0) for depth in depth_thresholds)
        self.rates_within_threshs = None    # defaultdict(float)

        self.missed_by_db = dict()
        self.var_num = None

    # def get_chrom_num(self):
    #     digits = [c for c in self.chrom if c.isdigit()]
    #     if digits:
    #         return int(''.join(digits))
    #     if 'M' in self.chrom:
    #         return 0
    #     if 'X' in self.chrom:
    #         return 23
    #     if 'Y' in self.chrom:
    #         return 24
    #     else:
    #         return 25
    #
    # @staticmethod
    # def get_order_key(r):
    #     return r.get_chrom_num(), r.get_start(), r.get_end(), r.gene_name

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_avg_depth(self):
        return self.avg_depth

    def get_std_dev(self):
        return self.std_dev

    def get_bases_by_depth(self):
        return self.bases_by_depth

    def get_size(self):
        if self.size is not None:
            return self.size
        if self.start and self.end:
            return abs(self.end - self.start)
        return None

    def add_bases_for_depth(self, depth, bases):
        if depth in self.bases_by_depth:
            err('Duplicated depth ' + str(depth) + ' for the region ' + str(self))
        else:
            self.bases_by_depth[depth] += bases

    def __str__(self):
        ts = [self.sample_name, self.chrom, self.start, self.end,
              self.gene_name, self.feature]
        return '"' + '\t'.join(map(str, ts)) + '"'

    def to_bed_str(self):
        return '\t'.join([self.chrom, str(self.start), str(self.end), self.gene_name] + '\n')

    def __repr__(self):
        return self.__str__()

    def calc_bases_within_threshs(self, depth_thresholds):
        if self.bases_within_threshs is not None:
            return self.bases_within_threshs

        if self.bases_by_depth is None:
            err('Error: self.bases_by_depth is None for ' + str(self))

        self.bases_within_threshs, self.rates_within_threshs = calc_bases_within_threshs(
            self.bases_by_depth, self.get_size(), depth_thresholds)

        return self.bases_within_threshs

    def calc_avg_depth(self):
        if self.avg_depth is not None:
            return self.avg_depth

        if self.bases_by_depth:
            depth_sum = sum(
                depth * bases
                for depth, bases
                in self.bases_by_depth.items())
            self.avg_depth = float(depth_sum) / self.get_size() if self.get_size() else None
            return self.avg_depth

    def calc_std_dev(self, avg_depth):
        if avg_depth is None:
            return None

        if self.std_dev is not None:
            return self.std_dev

        if self.bases_by_depth:
            sum_of_sq_var = sum(
                (depth - avg_depth) ** 2 * bases
                for depth, bases
                in self.bases_by_depth.items())
            sz = self.get_size()
            if sz and sz > 0:
                d = float(sum_of_sq_var) / float(sz)
                try:
                    self.std_dev = math.sqrt(d)
                except ValueError, e:
                    print 'float(sum_of_sq_var) =', float(sum_of_sq_var)
                    print 'float(sz) =', float(sz)
                    print 'd =', d
                    print self.sample_name, self.gene_name, self.chrom, ':', self.start, '-', self.end
                    # print 'math.sqrt(d) =', math.sqrt(d)
                    critical(str(e))
            return self.std_dev

    def calc_rate_within_normal(self, avg_depth):
        if avg_depth is None:
            return None

        if self.rate_within_normal is not None:
            return self.rate_within_normal

        return calc_rate_within_normal(self.bases_by_depth, avg_depth, self.get_size())

    def sum_up(self, depth_thresholds):
        self.calc_avg_depth()
        self.calc_std_dev(self.avg_depth)
        self.calc_bases_within_threshs(depth_thresholds)
        self.calc_rate_within_normal(self.avg_depth)
        return self.bases_within_threshs,\
               self.bases_by_depth, \
               self.avg_depth, \
               self.std_dev, \
               self.rate_within_normal

    def intersect(self, reg2):
        return self.chrom == reg2.chrom and \
               (reg2.start < self.start < reg2.end or
                self.start < reg2.start < self.end)


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
    def __init__(self, sample_name, gene_name, chrom=None, strand=None, feature='Gene-Exon', exon_num=None):
        Region.__init__(self,
            sample_name=sample_name, gene_name=gene_name, exon_num=exon_num, strand=strand,
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


def build_gene_objects_list(sample_name, features_bed, gene_keys_list):
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

    info('Building the Gene objects list based on target')
    if gene_keys_list:
        info()
        info('Init the Gene object dict')
        for gn, chrom in gene_keys_list:
            gene_by_name_and_chrom[(gn, chrom)] = GeneInfo(sample_name=sample_name, gene_name=gn, chrom=chrom)
        info('Processed ' + str(len(gene_keys_list)) + ' gene records -> ' + str(len(gene_by_name_and_chrom)) + ' uniq gene sybmols')

    info('Building the Gene objects list based on target')
    if features_bed and gene_by_name_and_chrom:
        # info()
        # info('Filtering exon bed file to have only gene records...')
        # exons_only_genes_bed = intermediate_fname(cnf, exons_bed, 'only_genes')
        # call(cnf, 'grep -w Gene ' + exons_bed, output_fpath=exons_only_genes_bed)
        # info('Saved genes to ' + exons_only_genes_bed)

        info()
        info('Setting start and end for the genes (based only on the target gene names found in the features list)')
        i = 0
        with open(features_bed) as f:
            for l in f:
                l = l.strip()
                if l and not l.startswith('#'):
                    chrom, start, end, symbol, _, strand, feature, biotype, transcript = l.split('\t')
                    if feature == 'Transcript':
                        gene_by_name_and_chrom[(symbol, chrom)].chrom = chrom
                        gene_by_name_and_chrom[(symbol, chrom)].start = int(start)
                        gene_by_name_and_chrom[(symbol, chrom)].end = int(end)
                        gene_by_name_and_chrom[(symbol, chrom)].size = int(end) - int(start)
                        gene_by_name_and_chrom[(symbol, chrom)].strand = strand
                        gene_by_name_and_chrom[(symbol, chrom)].biotype = biotype
                        gene_by_name_and_chrom[(symbol, chrom)].transcript_id = transcript
                        i += 1
        info('Processed ' + str(i) + ' genes')
        info()

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


def save_regions_to_bed(cnf, regions, bed_fpath, save_original_fields=False):
    if isfile(bed_fpath):
        if cnf.reuse_intermediate:
            verify_file(bed_fpath, is_critical=True)
            return bed_fpath
        else:
            os.remove(bed_fpath)

    with file_transaction(cnf.work_dir, bed_fpath) as tx_fpath:
        save_regions_to_bed_nocnf(regions, tx_fpath, save_original_fields)
    return bed_fpath


def save_regions_to_bed_nocnf(regions, bed_fpath, save_original_fields=False):
    info('Writing regions to ' + bed_fpath)

    with open(bed_fpath, 'w') as f:
        for r in regions:
            ts = [r.chrom, str(r.get_start()), str(r.get_end()), r.gene_name or '.']
            if save_original_fields:
                ts.extend(r.extra_fields)
            else:
                ts.append(r.feature or '.')
            f.write('\t'.join(ts) + '\n')

                # r.size, r.avg_depth, r.std_dev, r.percent_within_normal
            # f.write('\t'.join([
            #     region.chrom,
            #     str(region.start), str(region.end), str(region.avg_depth)]) + '\n')

    info('Saved to ' + bed_fpath)