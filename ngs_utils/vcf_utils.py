import re
from logging import critical
from cyvcf2 import VCF
from ngs_utils.file_utils import open_gzipsafe


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
    else:
        guessed_t_name, guessed_n_name = guess_sample_names(vcf_path)
        if not guessed_t_name:
            critical(f'Can\'t guess sample names from the VCF {vcf_path}')
        t_name = guessed_t_name
        n_name = guessed_n_name


    t_id = vcf_samples.index(t_name)
    if len(vcf_samples) >= 2:
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

