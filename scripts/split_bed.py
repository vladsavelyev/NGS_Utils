#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

from os.path import abspath, dirname, realpath, join
import sys
from source.utils import human_sorted


""" Input: Sorted and annotated BED file (i.e. at least 4 columns), File with list of key genes
    Output:
        BED file with no overlapped version of regions from input. In case of overlap scripts prefers 1) key gene,
        2) most upstream gene.

    Usage: python split_bed Input_BED_file [Key_genes_file] > Split_BED_file
"""


def _read_args(args):
    if len(args) < 1:
        log('Usage:\n')
        log('  ' + __file__ + ' Input_BED_file [Key_genes_file] > Split_BED_file\n')
        sys.exit(1)

    input_bed_fpath = abspath(args[0])
    log('Input: ' + input_bed_fpath)

    key_genes_fpath = '/ngs/reference_data/genomes/Hsapiens/common/az_key_genes.txt'
    if len(args) > 1:
        if args[1]:
            key_genes_fpath = abspath(args[1])
            log('Over-set key genes fpath: ' + key_genes_fpath)
    log()

    return input_bed_fpath, key_genes_fpath


def log(msg=''):
    sys.stderr.write(msg + '\n')


class Region:
    def __init__(self, entries):  # entries = line.strip().split('\t')
        self.chrom = entries[0]
        self.start = int(entries[1])
        self.end = int(entries[2])
        self.symbol = entries[3]
        self.rest = '\t'.join(entries[4:])

    def subregion(self, start=None, end=None):
        return Region([self.chrom, start if start else self.start, end if end else self.end, self.symbol, self.rest])

    def __str__(self):
        return '\t'.join([self.chrom, str(self.start), str(self.end), self.symbol, self.rest]) + '\n'

    def overlaps(self, other):  # assume 0-based coordinates!
        if self.chrom == other.chrom and self.end > other.start:
            return True
        return False

    def is_empty(self):
        return self.start == self.end

    @staticmethod
    def non_empty_list(regions):
        return [r for r in regions if not r.is_empty()]

    def split(self, other, key_genes):
        if self.end > other.end:  # other is fully included inside self
            if other.symbol in key_genes:  # note: two key genes can't overlap
                return Region.non_empty_list([self.subregion(end=other.start), other, self.subregion(start=other.end)])
            else:
                return Region.non_empty_list([self])
        else:
            if other.symbol in key_genes:  # note: two key genes can't overlap
                return Region.non_empty_list([self.subregion(end=other.start), other])
            else:
                return Region.non_empty_list([self, other.subregion(start=self.end)])

    def get_key(self):
        return '\t'.join([self.chrom, str(self.start), str(self.end)])

    def __lt__(self, other):
        # special case: chrM goes to the end
        if self.chrom != other.chrom and (self.chrom == 'chrM' or other.chrom == 'chrM'):
            return True if other.chrom == 'chrM' else False
        sorted_pair = human_sorted([self.get_key(), other.get_key()])
        if sorted_pair[0] == self.get_key() and sorted_pair[1] != self.get_key():
            return True
        return False


def main():
    input_bed_fpath, key_genes_fpath = _read_args(sys.argv[1:])

    key_genes = []
    with open(key_genes_fpath, 'r') as f:
        for line in f:
            key_genes.append(line.strip())

    input_regions = []
    with open(input_bed_fpath, 'r') as f:
        for line in f:
            if line.startswith('track') or line.startswith('browser'):
                continue
            if len(line.strip().split('\t')) < 4:
                log("ERROR: Input_BED_file should be annotated (i.e. has at least 4 columns)!")
                sys.exit(1)
            input_regions.append(Region(line.strip().split('\t')))

    prev_regions = []
    for cur_region in sorted(input_regions):
        # remove all previous regions before (i.e. not overlapping) cur_region
        while prev_regions:
            if not prev_regions[0].overlaps(cur_region):
                region = prev_regions.pop(0)
                sys.stdout.write(str(region))
            else:
                break
        if prev_regions:
            # splitting overlaps of cur_region with rest of prev_regions
            i = 0
            while i < len(prev_regions):
                region = prev_regions[i]
                if region.overlaps(cur_region):
                    split_regions = region.split(cur_region, key_genes)
                    if len(prev_regions) > i + 1 and split_regions[-1].overlaps(prev_regions[i + 1]):
                        prev_regions = prev_regions[:i] + split_regions[:-1] + prev_regions[i + 1:]
                        i += len(split_regions) - 1
                        cur_region = split_regions[-1]
                    else:  # no overlaps with rest of prev_regions
                        prev_regions = prev_regions[:i] + split_regions + prev_regions[i + 1:]
                        break
                else:
                    break
        else:
            prev_regions = [cur_region]
    for region in prev_regions:
        sys.stdout.write(str(region))


if __name__ == '__main__':
    main()
