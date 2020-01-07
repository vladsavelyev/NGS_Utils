from os.path import join, isfile

from natsort import natsort_keygen

from ngs_utils.file_utils import adjust_path, verify_file
from ngs_utils.logger import info, warn


class BaseSample:
    natsort_key = natsort_keygen()

    def __init__(self, name=None, dirpath=None, work_dir=None, bam=None, bed=None, vcf=None, genome=None,
                 targqc_dirpath=None, clinical_report_dirpath=None,
                 normal_match=None, sv_fpath=None, sv_bed=None,
                 l_fpath=None, r_fpath=None, phenotype=None, parent_project=None, **kwargs):
        self.name = name
        self.dirpath = dirpath
        self.work_dir = work_dir
        self.bam = bam
        self.counts_file = None
        self.l_fpath = l_fpath
        self.r_fpath = r_fpath
        self.is_wgs = False
        self.vcf = vcf
        self.phenotype = phenotype
        self.gender = None
        self.genome = genome
        self.var_dirpath = None
        self.normal_match = normal_match
        self.min_af = None
        self.sv_fpath = sv_fpath
        self.targqc_dirpath = targqc_dirpath
        self.clinical_html = None

        self.parent_project = parent_project
        self.is_rnaseq = None
        self.is_wgs = None

        self.min_allele_fraction = None
        self.coverage_interval = None
        self.variant_regions_bed = None
        self.coverage_bed = None

        self.batch = None
        self.batch_names = []
        self.phenotype = None

        for k, v in kwargs.items():
            self.__dict__[k] = v

    def __lt__(self, other):
        return self.key_to_sort() < other.key_to_sort()

    def key_to_sort(self):
        return BaseSample.natsort_key(self.name)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    @classmethod
    def load(cls, data):
        sample = cls(**data)
        sample.__dict__ = data
        return sample


class Batch:
    def __init__(self, name=None, normal=None, tumor=None, parent_project=None):
        self.name = name
        self.normal = normal
        self.tumor = tumor
        self.parent_project = parent_project

        self.somatic_vcf = None
        self.germline_vcf = None

    def is_paired(self):
        return self.normal and self.tumor

    def is_germline(self):
        return self.tumor.phenotype == 'germline'

    def __str__(self):
        return self.name

    def find_somatic_vcf(self):
        pass

    def find_germline_vcf(self):
        pass


class BcbioBatch(Batch):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def find_somatic_vcf(self, silent=False, caller=None):
        caller = caller or self.parent_project.somatic_caller

        # in datestamp. cwl-bcbio writes there
        vcf_cwl_fpath_gz = adjust_path(join(self.parent_project.date_dir, self.name + '-' + caller + '.vcf.gz'))
        # in datestamp. bcbio before 1.1.6
        vcf_old_fpath_gz = adjust_path(join(self.parent_project.date_dir, self.name + '-' + caller + '-annotated.vcf.gz'))
        # in sample dir. starting from bcbio 1.1.6, ~ Dec 2019
        vcf_fpath_gz = adjust_path(join(self.tumor.dirname, self.tumor.name + '-' + caller + '.vcf.gz'))

        if isfile(vcf_fpath_gz):
            verify_file(vcf_fpath_gz, is_critical=True)
            if not silent: info(f'Found somatic VCF in <final-dir>/<tumor-name>/<tumor-name>-{caller}.vcf.gz (conventional bcbio): ' + vcf_fpath_gz)
            self.somatic_vcf = vcf_fpath_gz

        elif isfile(vcf_old_fpath_gz):
            verify_file(vcf_old_fpath_gz, is_critical=True)
            if not silent: info(f'Found somatic VCF in <date-dir>/<batch>-{caller}-annotated.vcf.gz (bcbio < v1.1.6)): ' + vcf_old_fpath_gz)
            self.somatic_vcf = vcf_old_fpath_gz

        elif isfile(vcf_cwl_fpath_gz):
            verify_file(vcf_cwl_fpath_gz, is_critical=True)
            if not silent: info(f'Found somatic VCF in project/<batch>-{caller}.vcf.gz (CWL bcbio): ' + vcf_cwl_fpath_gz)
            self.somatic_vcf = vcf_cwl_fpath_gz

        elif not silent:
            warn(f'Could not find somatic variants files for batch {self.name}, caller {caller} neither as '
                 f'{self.parent_project.final_dir}/<tumor-name>/<tumor-name>-{caller}.vcf.gz (conventional bcbio), nor as '
                 f'{self.parent_project.date_dir}/<batch>-{caller}-annotated.vcf.gz (bcbio < v1.1.6), nor as '
                 f'project/<batch>-{caller}.vcf.gz (CWL bcbio).')

    def find_germline_vcf(self, silent=False, caller=None):
        caller = caller or self.parent_project.germline_caller

        # in sample dir. starting from bcbio 1.1.6, ~ Dec 2019
        vcf_fpath_gz = adjust_path(join(self.parent_project.date_dir, f'{self.normal.name}-germline-{caller}.vcf.gz'))
        # in datestamp. bcbio before 1.1.6
        vcf_old_fpath_gz = adjust_path(join(self.parent_project.date_dir, f'{self.normal.name}-germline-{caller}-annotated.vcf.gz'))

        if isfile(vcf_fpath_gz):
            verify_file(vcf_fpath_gz, is_critical=True)
            if not silent: info(f'Found germline VCF in <date-dir>/<normal-name>-germline-{caller}.vcf.gz: ' + vcf_fpath_gz)
            self.germline_vcf = vcf_fpath_gz

        elif isfile(vcf_old_fpath_gz):
            verify_file(vcf_old_fpath_gz, is_critical=True)
            if not silent: info(f'Found germline VCF in <date-dir>/<normal-name>-germline-{caller}-annotated.vcf.gz (bcbio < v1.1.6)): ' + vcf_old_fpath_gz)
            self.germline_vcf = vcf_old_fpath_gz

        elif not silent:
            warn(f'Could not find germline variants files for batch {self.name}, caller {caller} neither as '
                 f'<date-dir>/<normal-name>-germline-{caller}.vcf.gz, nor as '
                 f'<date-dir>/<normal-name>-germline-{caller}-annotated.vcf.gz (bcbio < v1.1.6)')


# class Caller:
#     def __init__(self, name=None, is_germline=False):
#         self.name = name
#         self.is_germline = is_germline
#         self.samples = []
#
#     def __str__(self):
#         return self.name