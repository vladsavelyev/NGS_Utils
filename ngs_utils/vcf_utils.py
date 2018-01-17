import re
from ngs_utils.file_utils import open_gzipsafe


def get_sample_names(vcf_path):
    return get_sample_ids(vcf_path)


def get_tumor_sample_name(vcf_path):
    return get_sample_ids(vcf_path, return_names=True)[0]


def get_normal_sample_name(vcf_path):
    return get_sample_ids(vcf_path, return_names=True)[1]


def get_tumor_sample_id(vcf_path):
    return get_sample_ids(vcf_path)[0]


def get_normal_sample_id(vcf_path):
    return get_sample_ids(vcf_path)[1]


def get_sample_ids(vcf_path, return_names=False):
    """ Finds tumor and control sample names/ids from a bcbio-derived VCF,
        and returns a tuple (tumor, control)

        1. If ##SAMPLE header field is found, use that for the name.
        2. Otherwise, return the first (for tumor) and second (for normal)
           sample in #CHROM header.
    """
    tumor_name, control_name = None, None
    tumor_id, control_id = None, None

    with open_gzipsafe(vcf_path) as f:
        for line in f:
            m = re.match(r'^##SAMPLE=<ID=(?P<name>\S+),Genomes=Tumor>$', line)
            if m:
                tumor_name = m.group('name')
            m = re.match(r'^##SAMPLE=<ID=(?P<name>\S+),Genomes=Germline>$', line)
            if m:
                control_name = m.group('name')

            if tumor_name and return_names:
                return tumor_name, control_name

            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                if tumor_name:
                    tumor_id = samples.index(tumor_name)
                    if control_name:
                        control_id = samples.index(control_name)
                else:
                    tumor_id, control_id = 0, None
                    tumor_name = samples[tumor_id]
                    if len(samples) > 1:
                        control_id = 1
                        control_name = samples[control_id]

                if return_names:
                    return tumor_name, control_name
                else:
                    return tumor_id, control_id
    raise ValueError


