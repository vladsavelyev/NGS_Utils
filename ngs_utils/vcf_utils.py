import re
import subprocess
import sys
from logging import critical
from ngs_utils.file_utils import open_gzipsafe
from cyvcf2 import VCF, Writer
from ngs_utils.call_process import run_simple
from ngs_utils.file_utils import get_ungz_gz


def get_sample_names(
        vcf_path,
        provided_tumor_name=None,
        provided_normal_name=None):
    return get_sample_ids(
        vcf_path,
        provided_t_name=provided_tumor_name,
        provided_n_name=provided_normal_name,
        return_names=True
    )


def get_tumor_sample_name(
        vcf_path,
        provided_tumor_name=None,
        provided_normal_name=None):
    return get_sample_ids(
        vcf_path,
        provided_t_name=provided_tumor_name,
        provided_n_name=provided_normal_name,
        return_names=True)[0]


def get_normal_sample_name(
        vcf_path,
        provided_tumor_name=None,
        provided_normal_name=None):
    return get_sample_ids(
        vcf_path,
        provided_t_name=provided_tumor_name,
        provided_n_name=provided_normal_name,
        return_names=True)[1]


def get_tumor_sample_id(
        vcf_path,
        provided_tumor_name=None,
        provided_normal_name=None):
    return get_sample_ids(
        vcf_path,
        provided_t_name=provided_tumor_name,
        provided_n_name=provided_normal_name)[0]


def get_normal_sample_id(
        vcf_path,
        provided_tumor_name=None,
        provided_normal_name=None):
    return get_sample_ids(
        vcf_path,
        provided_t_name=provided_tumor_name,
        provided_n_name=provided_normal_name)[1]


def get_sample_ids(
        vcf_path,
        provided_t_name=None,
        provided_n_name=None,
        return_names=False):

    t_id, n_id = None, None
    t_name, n_name = None, None

    vcf_samples = VCF(vcf_path).samples

    if provided_t_name:
        assert provided_t_name in vcf_samples,\
            f'Tumor sample name {provided_t_name} is not in VCF {vcf_path}. Found: {vcf_samples}'
        t_name = provided_t_name
    if provided_n_name:
        assert provided_n_name in vcf_samples,\
            f'Normal sample name {provided_n_name} is not in VCF {vcf_path}. Found: {vcf_samples}'
        n_name = provided_n_name
    guessed_t_name, guessed_n_name = guess_sample_names(vcf_path)
    if not t_name:
        if not guessed_t_name:
            critical(f'Can\'t guess tumor sample name from the VCF {vcf_path}')
        t_name = guessed_t_name
    if not n_name:
        if not guessed_n_name:
            critical(f'Can\'t guess normal sample name from the VCF {vcf_path}')
        n_name = guessed_n_name

    t_id = vcf_samples.index(t_name)
    if n_name and len(vcf_samples) >= 2:
        n_id = vcf_samples.index(n_name)

    if return_names:
        return t_name, n_name
    else:
        return t_id, n_id


def guess_sample_names(vcf_path):
    """ Finds tumor and control sample names from a bcbio-derived VCF,
        and returns a tuple (tumor, control)

        1. If ##SAMPLE header field is found, use that for the names
        2. If ##DRAGENCommandLine is found, with --RGSM parameters, use them
        3. Otherwise, we don't guess and rathe rely on setting them explicitly
           (see -tn and -nn options in scripts and provided_tumor_name and
           provided_normal_name in get_sample_ids).

        Returns: tumor_name, normal_name
    """
    t_name, n_name = None, None

    with open_gzipsafe(vcf_path) as f:
        for line in f:
            # bcbio?
            m = re.match(r'^##PEDIGREE=<Derived=(?P<t_name>\S+),Original=(?P<n_name>\S+)>$', line)
            if m:
                t_name = m.group('t_name')
                n_name = m.group('n_name')
            else:
                m = re.match(r'^##SAMPLE=<ID=(?P<name>\S+),Genomes=Tumor>$', line)
                if m:
                    t_name = m.group('name')
                m = re.match(r'^##SAMPLE=<ID=(?P<name>\S+),Genomes=Germline>$', line)
                if m:
                    n_name = m.group('name')

            # dragen?
            m = re.match(r'^##DRAGENCommandLine=<ID=dragen,.*'
                         r'--RGSM\s+(?P<n_name>\S+)\s+.*'
                         r'--RGSM-tumor\s+(?P<t_name>\S+)\s+.*', line)
            if m:
                t_name = m.group('t_name')
                n_name = m.group('n_name')
            else: # dragen single-sample?
                m = re.match(r'^##DRAGENCommandLine=<ID=dragen,.*'
                             r'--RGSM\s+(?P<s_name>\S+)\s+.*', line)
                if m:
                    t_name = m.group('s_name')

    return t_name, n_name


def add_cyvcf2_hdr(vcf, id, number, type, descr, new_header=None, hdr='INFO'):
    if new_header:
        new_header.append(f'##{hdr}=<ID={id},Number={number},Type={type},Description="{descr}">')
    if hdr == 'INFO':
        vcf.add_info_to_header({'ID': id, 'Type': type, 'Number': number, 'Description': descr})
    elif hdr == 'FORMAT':
        vcf.add_format_to_header({'ID': id, 'Type': type, 'Number': number, 'Description': descr})
    else:
        critical(f'Unknown VCF header: {hdr}. Supported: INFO, FORMAT')


def count_vars(vcf_path, filter_col=None, bcftools_filter_expr=None):
    cmd = f'cat {vcf_path} | '
    if bcftools_filter_expr:
        cmd += f'bcftools filter {bcftools_filter_expr} |'
    if filter_col:
        cmd += f'bcftools view -f {filter_col} | '
    cmd += 'bcftools view -H | wc -l'
    return int(subprocess.check_output(cmd, shell=True).strip())


def vcf_contains_field(vcf_path, field, col=None):
    # col is FILTER, FORMAT, INFO, or None (=any of three)
    if col is None and '/' in field:
        col, field = field.split('/')
    if col is not None:
        return f'##{col}=<ID={field},' in VCF(vcf_path).raw_header
    return f'=<ID={field},' in VCF(vcf_path).raw_header


def iter_vcf(input_file, output_file, proc_rec, proc_hdr=None, postproc_hdr=None, **kwargs):
    """
    :param input_file: path to input VCF file
    :param output_file: path to output VCF file (can be .vcf or .vcf.gz, but it will always bgzip/tabix and write with .vcf.gz extention)
    :param proc_rec: a function to process a single cyvcf Record object. Returns either a (new) Record object to write, or None to indicate that the record should be discarded
    :param proc_hdr: a function to process cyvcf object once (i.e. to add values to the header with vcf.add_info_to_header, etc)
    :param postproc_hdr: a function to postprocess finalized header string (vcf.rawheader), e.g. in order to remove values
    :param kwargs: any paramters to pass directly into proc_rec
    """
    vcf = VCF(input_file, gts012=True)
    if proc_hdr is not None:
        proc_hdr(vcf)

    # w = None
    if output_file is not None:
        out_ungz, out_gz = get_ungz_gz(output_file)
        # w = Writer(out_ungz, vcf)
        # w.write_header()
        w = open(out_ungz, 'w')
    else:
        # sys.stdout.write(vcf.raw_header)
        w = sys.stdout

    header = vcf.raw_header
    if postproc_hdr is not None:
        header = postproc_hdr(header)
    w.write(header)

    for rec in vcf:
        if proc_rec:
            rec_res = proc_rec(rec, vcf, **kwargs)
            if rec_res is not None:
                # if w is not None:
                #     sys.stderr.write('Writing record', rec_res, '\n')
                #     w.write_record(rec_res)
                # else:
                #     print(rec_res)
                # sys.stderr.write(f'Writing record {rec_res}\n')
                w.write(f'{rec_res}')

    sys.stderr.write(f'Finished writing {output_file}\n')
    vcf.close()
    if output_file is not None:
        w.close()
        run_simple(f'bgzip -f {out_ungz} && tabix -f -p vcf {out_gz}')
        sys.stderr.write(f'Compressed {output_file}\n')


def iter_vcf__pysam(input_file, proc_rec=None, proc_hdr=None, output_file=None):
    import pysam
    import sys

    vcf = pysam.VariantFile(input_file)
    if output_file:
        w = open(output_file, 'w')
    else:
        w = sys.stdout

    # Header
    if proc_hdr is not None:
        proc_hdr(vcf)
    w.write(str(vcf.header))

    # Records
    for rec in vcf:
        if proc_rec:
            rec_res = proc_rec(rec)
            if rec_res is not None:
                print(rec_res)
                w.write(str(rec_res))

    vcf.close()

    if output_file:
        w.close()
        out_ungz, out_gz = get_ungz_gz(output_file)
        run_simple(f'bgzip -f {out_ungz} && tabix -f -p vcf {out_gz}')

