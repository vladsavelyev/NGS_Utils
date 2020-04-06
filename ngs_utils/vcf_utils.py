import re
from logging import critical

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

    guessed_t_name, guessed_n_name, guessed_t_id, guessed_n_id = \
        guess_sample_names_and_ids(vcf_path)
    sample_names = [guessed_t_name, guessed_n_name]
    sample_ids = [guessed_t_id, guessed_n_id]

    t_name, n_name = None, None

    if provided_t_name:
        t_name = provided_t_name
        # guessing normal if only provided tumor name:
        if not provided_n_name and len(sample_names) >= 2:
            n_name = [sn for sn in sample_names if sn != t_name][0]
    elif not t_name:
        t_name = guessed_t_name

    if provided_n_name:
        n_name = provided_n_name
        # guessing tumor if only provided tumor name:
        if not provided_t_name and len(sample_names) >= 2:
            t_name = [sn for sn in sample_names if sn != n_name][0]
    elif not n_name:
        n_name = guessed_n_name

    n_id = None
    if n_name:
        assert n_name in sample_names, 'n_name: ' + str(n_name) + ', sample_names: ' + str(sample_names)
        n_id = sample_ids[sample_names.index(n_name)]
    assert t_name in sample_names, 't_name: ' + str(t_name) + ', sample_names: ' + str(sample_names)
    t_id = sample_ids[sample_names.index(t_name)]

    if return_names:
        return t_name, n_name
    else:
        return t_id, n_id


def guess_sample_names_and_ids(vcf_path):
    """ Finds tumor and control sample names/ids from a bcbio-derived VCF,
        and returns a tuple (tumor, control)

        1. If ##SAMPLE header field is found, use that for the name.
        2. Otherwise, return the first (for tumor) and second (for normal)
           sample in #CHROM header. It would usually be wrong, so better rely
           on setting them explicitly (see -tn and -nn options and provided_tumor_name
           and provided_normal_name in get_sample_ids).

        Returns: tumor_name, normal_name, tumor_id, normal_id
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

            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                tumor_id, normal_id = None, None
                if t_name:
                    tumor_id = samples.index(t_name)
                if n_name and n_name in samples:
                    normal_id = samples.index(n_name)

                if tumor_id is None and normal_id is not None:
                    tumor_id = 1 if normal_id == 0 else 0

                elif normal_id is None and tumor_id is not None and len(samples) > 1:
                    normal_id = 1 if tumor_id == 0 else 0

                elif normal_id is None and tumor_id is None:
                    if len(samples) == 1:
                        tumor_id = 0
                    else:
                        raise ValueError(f'Can\'t guess sample names from the VCF {vcf_path}')

                if t_name is None:
                    t_name = samples[tumor_id]
                if n_name is None and len(samples) > 1:
                    n_name = samples[normal_id]

                return t_name, n_name, tumor_id, normal_id
    raise ValueError(f'Can\'t guess sample names from the VCF {vcf_path}')


def add_cyvcf2_hdr(vcf, id, number, type, descr, new_header=None, hdr='INFO'):
    if new_header:
        new_header.append(f'##{hdr}=<ID={id},Number={number},Type={type},Description="{descr}">')
    if hdr == 'INFO':
        vcf.add_info_to_header({'ID': id, 'Type': type, 'Number': number, 'Description': descr})
    elif hdr == 'FORMAT':
        vcf.add_format_to_header({'ID': id, 'Type': type, 'Number': number, 'Description': descr})
    else:
        critical(f'Unknown VCF header: {hdr}. Supported: INFO, FORMAT')

