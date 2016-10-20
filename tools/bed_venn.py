
import sys
from optparse import OptionParser
import subprocess
from os.path import dirname, basename, join, splitext, isfile
import json
import itertools

from Utils.bed_utils import verify_bed
from Utils.file_utils import adjust_path, safe_mkdir
from Utils.logger import critical

# beds_dirpath = '/Users/vlad/googledrive/az/reference_data/genomes/Hsapiens/hg19/bed/intersection/merged'
# work_dirpath = '/Users/vlad/googledrive/az/reference_data/genomes/Hsapiens/hg19/bed/intersection/merged/intersections'


def call(cmdl):
    print cmdl
    subprocess.call(cmdl, shell=True)


def check_output(cmdl):
    print cmdl
    return subprocess.check_output(cmdl, shell=True)


def bedsize(bed):
    size = check_output("cat " + bed + " | awk -F'\\t' 'BEGIN{ SUM=0 }{ SUM+=$3-$2 }END{ print SUM }'")
    print 'Size of ' + basename(bed) + ': ' + size
    return int(size)


total_calls = 0


def intersect_pair(work_dirpath, bed1, bed2):
    bed1, bed2 = sorted([bed1, bed2])
    output_fpath = join(work_dirpath, splitext(basename(bed1))[0] + '__' + basename(bed2))
    if not isfile(output_fpath):
        print 'intersect_pair: ' + splitext(basename(bed1))[0] + ' and ' + splitext(basename(bed2))[0]
        call('/usr/local/bin/bedtools intersect -a ' + bed1 + ' -b ' + bed2 + ' > ' + output_fpath)
        # global total_calls
        # total_calls += 1
    return output_fpath


def calc_set_intersection(work_dirpath, sorted_beds_subset, intersection_bed_by_subset):
    if len(sorted_beds_subset) == 1:
        return sorted_beds_subset[0]
    if tuple(sorted_beds_subset) in intersection_bed_by_subset:
        return intersection_bed_by_subset[tuple(sorted_beds_subset)]

    intersection_bed = None
    for bed in sorted_beds_subset:
        remaining_subset = [b for b in sorted_beds_subset if b != bed]  # sorted_beds_subset.index(b) > sorted_beds_subset.index(bed)]
        if not remaining_subset:
            continue
        print 'comparing ' + basename(bed) + ' and ' + str([basename(k) for k in remaining_subset])
        subset_intersection_bed = intersection_bed_by_subset.get(tuple(sorted(remaining_subset)))
        if not subset_intersection_bed:
            subset_intersection_bed = calc_set_intersection(work_dirpath, remaining_subset, intersection_bed_by_subset)
        intersection_bed = intersect_pair(work_dirpath, subset_intersection_bed, bed)
        intersection_bed_by_subset[tuple(sorted(remaining_subset))] = subset_intersection_bed
    intersection_bed_by_subset[tuple(sorted_beds_subset)] = intersection_bed
    return intersection_bed


def save_venn_diagram_data(size_by_set, label_by_set):
    data = []
    for venn_set, size in size_by_set.iteritems():
        set_info = dict()
        set_info['size'] = size
        if isinstance(venn_set, tuple):
            set_info['sets'] = list(venn_set)
        else:
            set_info['sets'] = [venn_set]
        if isinstance(venn_set, int):
            set_info['label'] = label_by_set[venn_set]
        data.append(set_info)
    return json.dumps(sorted(data, key=lambda x: x['sets']))


def main(bed_fpaths, output_dir):
    # beds = tuple([join(beds_dirpath, bed) for bed in [
    #     'AZ.bed',
    #     # 'CRE.bed',
    #     # 'IDT.bed',
    #     # 'Med.bed',
    #     'One.bed',
    #     'One.hg38tohg19.bed',
    #     # 'V6.bed',
    # ]])

    safe_mkdir(output_dir)
    work_dirpath = safe_mkdir(join(output_dir))

    intersection_bed_by_subset = dict()
    calc_set_intersection(work_dirpath, sorted(bed_fpaths), intersection_bed_by_subset)

    intersection_size_by_subset = dict()
    label_by_subset = dict()
    for bed_set, intersection_bed in intersection_bed_by_subset.items():
        bed_set = tuple([splitext(basename(b))[0] for b in bed_set])
        intersection_size_by_subset[bed_set] = bedsize(intersection_bed)
        label_by_subset[bed_set] = basename(splitext(intersection_bed)[0])
        print (str(bed_set) + ': ' + basename(intersection_bed) + ', size: ' + str(intersection_size_by_subset[bed_set]))

    json_txt = save_venn_diagram_data(intersection_size_by_subset, label_by_subset)
    print json_txt

    # TODO: generate HTML


if __name__ == '__main__':
    parser = OptionParser(usage='Usage: ' + basename(__file__) + ' bed1 bed2 ... -o results_dir')
    parser.add_option('-o', '--output-dir', dest='output_dir')
    (opts, args) = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        critical(parser.usage)
    if not opts.output_fpath:
        critical(parser.usage)

    bed_fpaths = [verify_bed(bed) for bed in args]
    output_dir = adjust_path(opts.output_dir)

    main(bed_fpaths, output_dir)


# sets = []
# for i in range(2, len(beds) + 1):
#     sets.append(list(itertools.combinations(beds, i)))


