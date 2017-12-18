import re
from ngs_utils.file_utils import open_gzipsafe


def get_tumor_sample_name(vcf_path):
    return _get_sample_id(vcf_path, 'Tumor', return_name=True)


def get_normal_sample_name(vcf_path):
    return _get_sample_id(vcf_path, 'Germline', return_name=True)


def get_tumor_sample_id(vcf_path):
    return _get_sample_id(vcf_path, 'Tumor')


def get_normal_sample_id(vcf_path):
    return _get_sample_id(vcf_path, 'Germline')


def _get_sample_id(vcf_path, key='Tumor', return_name=False):
    """ Finds tumor or control sample name/id from a bcbio-derived VCF.
        1. If ##SAMPLE header field is found, use that for the name.
        2. Otherwise, return the first (for tumor) or second (for normal)
        sample in #CHROM header.
    """
    assert key in ['Tumor', 'Germline']
    sample_name = None
    with open_gzipsafe(vcf_path) as f:
        for line in f:
            m = re.match(r'^##SAMPLE=<ID=(?P<name>\S+),Genomes=' + key + '>$', line)
            if m:
                sample_name = m.group('name')
                if return_name:
                    return sample_name

            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                if sample_name:
                    sample_id = samples.index(sample_name)
                else:
                    if key == 'Tumor':
                        sample_id = 0
                    else:
                        sample_id = 1

                if return_name:
                    return samples[sample_id]
                else:
                    return sample_id
    raise ValueError

