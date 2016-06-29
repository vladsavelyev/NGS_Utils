import os
from collections import OrderedDict
from os.path import isfile, join, abspath, basename, dirname, getctime, getmtime, splitext, realpath
from subprocess import check_output

from Utils import call_process
from Utils.file_utils import intermediate_fname, iterate_file, splitext_plus, verify_file, adjust_path, add_suffix, \
    safe_mkdir, file_transaction, which
from Utils.logger import info, critical, warn, err, debug
from Utils import reference_data as ref
from Utils.Region import SortableByChrom, get_chrom_order
from Utils.utils import md5


def markdup_sam(cnf, in_sam_fpath, samblaster=None):
    """Perform non-stream based deduplication of SAM input files using samblaster.
    """
    samblaster = samblaster or which('samblaster')
    if not samblaster:
        warn('No samblaster, can\'t mark duplicates.')
        return None

    out_sam_fpath = add_suffix(in_sam_fpath, 'markdup')
    tmp_fpath = join(cnf.work_dir, splitext_plus(basename(in_sam_fpath))[0] + '_markdup')
    safe_mkdir(dirname(tmp_fpath))
    cmdline = '{samblaster} -i {in_sam_fpath} -o {out_sam_fpath}'.format(**locals())
    res = call_process.run(cmdline, output_fpath=out_sam_fpath, stdout_to_outputfile=False)
    return out_sam_fpath


def markdup_bam(cnf, in_bam_fpath, bammarkduplicates=None):
    """Perform non-stream based deduplication of BAM input files using biobambam.
    """
    bammarkduplicates = bammarkduplicates or which('bammarkduplicates')
    if not bammarkduplicates:
        warn('No biobambam bammarkduplicates, can\'t mark duplicates.')
        return None

    out_bam_fpath = add_suffix(in_bam_fpath, 'markdup')
    tmp_fpath = join(cnf.work_dir, splitext_plus(basename(in_bam_fpath))[0] + '_markdup')
    safe_mkdir(dirname(tmp_fpath))
    cmdline = '{bammarkduplicates} tmpfile={tmp_fpath} I={in_bam_fpath} O={out_bam_fpath}'.format(**locals())
    res = call_process.run(cmdline, output_fpath=out_bam_fpath, stdout_to_outputfile=False)
    return out_bam_fpath


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


def remove_comments(work_dir, bed_fpath, reuse=False):
    def f(l, i):
        if not l.startswith('#'):
            return l
        else:
            return None
    return iterate_file(work_dir, bed_fpath, f, suffix='rmcmt', reuse=reuse)


def extract_gene_names_and_filter_exons(work_dir, target_bed, features_bed, reuse=False):
    gene_key_set = set()
    gene_key_list = []
    genes_not_in_refseq = set()

    debug()
    debug('Getting gene list')

    # if genes_fpath:
    #     with open(genes_fpath) as f:
    #         gene_key_list = [g.strip() for g in f.read().split('\n') if g]
    #         gene_key_set = set(gene_key_list)
    #     info('Using genes from ' + genes_fpath + ', filtering exons and amplicons with this genes.')
    #     if target_bed:
    #         target_bed = filter_bed_with_gene_set(cnf, target_bed, gene_key_set)
    #     if exons_bed:
    #         exons_bed = filter_bed_with_gene_set(cnf, exons_bed, gene_key_set)
    #         exons_no_genes_bed = filter_bed_with_gene_set(cnf, exons_no_genes_bed, gene_key_set)
    # else:

    if target_bed:
        debug()
        gene_key_set, gene_key_list = get_gene_keys(target_bed)
        debug('Using genes from the target ' + target_bed)
        if features_bed:
            debug('Trying filtering exons with these ' + str(len(gene_key_list)) + ' genes.')
            features_filt_bed, genes_in_refseq = filter_bed_with_gene_set(work_dir, features_bed, gene_key_set, suffix='target_genes_1st_round', reuse=reuse)
            if not genes_in_refseq:
                debug()
                warn('No gene symbols from the target BED file was found in the RefSeq features. Re-annotating target...')
                target_bed = annotate(target_bed, features_bed, add_suffix(target_bed, 'ann'), reuse=reuse)
                #info('Merging regions within genes...')
                #target_bed = group_and_merge_regions_by_gene(cnf, target_bed, keep_genes=False)
                # debug('Sorting amplicons_bed by (chrom, gene_name, start)')
                # target_bed = sort_bed(work_dir, target_bed)
                debug('Getting gene names again...')
                gene_key_set, gene_key_list = get_gene_keys(target_bed)
                debug()
                debug('Using genes from the new amplicons list, filtering features with this genes again.')
                features_filt_bed, genes_in_refseq = filter_bed_with_gene_set(work_dir, features_bed, gene_key_set, suffix='target_genes_2nd_round', reuse=reuse)
                if not genes_in_refseq:
                    critical('No gene symbols from the target BED file was found in the RefSeq features.')
            features_bed = features_filt_bed
            genes_not_in_refseq = gene_key_set - genes_in_refseq

    elif features_bed:
        info()
        info('No target (probably WGS), getting the gene names from the full features list...')
        gene_key_set, gene_key_list = get_gene_keys(features_bed)
        genes_not_in_refseq = set()
    info()

    return gene_key_set, gene_key_list, genes_not_in_refseq, target_bed, features_bed


def calc_region_number(bed_fpath):
    with open(bed_fpath) as f:
        return sum(1 for l in f if l.strip() and not l.strip().startswith('#'))


def get_gene_keys(bed_fpath, chrom_index=0, gene_index=3):
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

    def fn(l, i):
        if l:
            fs = l.split('\t')
            if len(fs) < 4:
                return None
            new_gns = []
            c = fs[0]
            for g in fs[3].split(','):
                if (g, c) in gene_keys_set:
                    new_gns.append(g)
                    met_genes.add((g, c))
            if new_gns:
                return l.replace(fs[3], ','.join(new_gns))

    res_fpath = iterate_file(work_dir, bed_fpath, fn, suffix=suffix or 'filt_genes', check_result=False)
    return res_fpath, met_genes


def sort_bed(input_bed_fpath, output_bed_fpath=None, work_dir=None, fai_fpath=None, genome=None, reuse=False):
    input_bed_fpath = verify_bed(input_bed_fpath)
    output_bed_fpath = adjust_path(output_bed_fpath) if output_bed_fpath \
        else intermediate_fname(work_dir, input_bed_fpath, 'sorted')

    class Region(SortableByChrom):
        def __init__(self, chrom, start, end, other_fields, chrom_ref_order):
            SortableByChrom.__init__(self, chrom, chrom_ref_order)
            self.start = start
            self.end = end
            self.chrom_ref_order = chrom_ref_order
            self.other_fields = tuple(other_fields)

        def get_key(self):
            return self.chrom_ref_order, self.start, self.end, self.other_fields

    regions = []

    if fai_fpath:
        fai_fpath = verify_file(fai_fpath)
    elif genome:
        fai_fpath = verify_file(ref.get_fai(genome))
    else:
        critical('fai or genome build name must be specified')
    chr_order = get_chrom_order(fai_fpath=fai_fpath)

    debug('Sorting regions in ' + input_bed_fpath)
    if reuse and isfile(output_bed_fpath) and verify_bed(output_bed_fpath):
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


def total_merge_bed(cnf, bed_fpath):
    bedops = 'bedops'
    if bedops:
        cmdline = '{bedops} --merge {bed_fpath}'.format(**locals())
        output_fpath = intermediate_fname(cnf, bed_fpath, 'total_merged')
        call_process.run(cmdline, output_fpath)
        return output_fpath
    else:
        bedtools = 'bedtools'
        cmdline = '{bedtools} merge -i {bed_fpath}'.format(**locals())
        output_fpath = intermediate_fname(cnf, bed_fpath, 'total_merged')
        call_process.run(cnf, cmdline, output_fpath)
        return output_fpath


def calc_sum_of_regions(bed_fpath):
    total_bed_size = 0

    with open(bed_fpath) as f:
        for l in f:
            l = l.strip()
            if not l.startswith('#'):
                start, end = [int(f) for f in l.split('\t')[1:3]]
                total_bed_size += end - start

    return total_bed_size


def get_total_bed_size(cnf, bed_fpath):
    merged_bed = total_merge_bed(cnf, bed_fpath)
    return calc_sum_of_regions(merged_bed)


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


def verify_bam(fpath, description='', is_critical=False, silent=False):
    if not verify_file(fpath, description, is_critical=is_critical, silent=silent):
        return None

    fpath = adjust_path(fpath)

    logfn = critical if is_critical else err
    if not fpath.endswith('.bam'):
        logfn('The file ' + fpath + ' is supposed to be BAM but does not have the .bam '
            'extension. Please, make sure you pass proper file.')
        return None

    textchars = ''.join(map(chr, [7, 8, 9, 10, 12, 13, 27] + range(0x20, 0x100)))
    is_binary_string = lambda baitiki: bool(baitiki.translate(None, textchars))
    if not is_binary_string(open(fpath).read(3)):
        logfn('The BAM file ' + fpath + ' must be a binary file.')
        return None

    return fpath


def verify_bed(fpath, description='', is_critical=False, silent=False):
    if not verify_file(fpath, description, is_critical=is_critical, silent=silent):
        return None

    fpath = adjust_path(fpath)

    error = BedFile(fpath).checkformat()
    if error:
        fn = critical if is_critical else err
        fn('Error: incorrect bed file format (' + fpath + '): ' + str(error) + '\n')
        return None

    return fpath


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

        fd = file(self.filename)

        line = fd.readline()
        fields = line.split('\t')
        lc = 1
        error = ''

        # Checks that the two columns on the right contain integer values
        try:
            # Parses each line and checks that there are at least 3 fields, the two on the right containing integer values and being the right one
            # greater than the left one
            while line <> '' and len(fields) > 2 and int(fields[1]) < int(fields[2]):
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