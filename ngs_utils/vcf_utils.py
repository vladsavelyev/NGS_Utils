import re
from ngs_utils.file_utils import open_gzipsafe


def get_sample_names(vcf_path):
    return get_sample_ids(vcf_path, return_names=True)


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
                tumor_id, control_id = None, None
                if tumor_name:
                    tumor_id = samples.index(tumor_name)
                if control_name and control_name in samples:
                    control_id = samples.index(control_name)

                if tumor_id is None and control_id is not None:
                    tumor_id = 1 if control_id == 0 else 0

                elif control_id is None and tumor_id is not None and len(samples) > 1:
                    control_id = 1 if tumor_id == 0 else 0

                elif control_id is None and tumor_id is None:
                    tumor_id = 0
                    if len(samples) > 1:
                        control_id = 1

                if tumor_name is None:
                    tumor_name = samples[tumor_id]
                    if control_name is None and len(samples) > 1:
                        control_name = samples[control_id]

                if return_names:
                    return tumor_name, control_name
                else:
                    return tumor_id, control_id
    raise ValueError




