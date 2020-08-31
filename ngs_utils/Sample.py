from natsort import natsort_keygen
from collections import Iterable


class BaseProject:
    def __init__(self, input_dir=None, **kwargs):
        self.dir = input_dir
        self.project_name = None
        self.samples = []
        self.somatic_caller = None
        self.germline_caller = None
        self.genome_build = None
        self.batch_by_name = dict()
        self.is_rnaseq = False
        self.is_wgs = True


class BaseSample:
    natsort_key = natsort_keygen()

    def __init__(self, name=None, dirpath=None, work_dir=None, bam=None, bed=None, genome=None,
                 targqc_dirpath=None, clinical_report_dirpath=None,
                 normal_match=None, sv_fpath=None, sv_bed=None,
                 l_fpath=None, r_fpath=None, phenotype=None, rgid=None,
                 parent_project=None, **kwargs):
        self.name = name
        self.rgid = rgid if rgid is not None else name
        self.dirpath = dirpath
        self.work_dir = work_dir
        self.bam = bam
        self.counts_file = None
        self.l_fpath = l_fpath
        self.r_fpath = r_fpath
        self.is_wgs = False
        self.phenotype = phenotype
        self.gender = None
        self.genome = genome
        self.var_dirpath = None
        self.normal_match = normal_match
        self.min_af = None
        self.sv_fpath = sv_fpath
        self.targqc_dirpath = targqc_dirpath
        self.clinical_html = None

        self.batches = []
        self.batch_names = []
        self.phenotype = None
        self.is_rnaseq = None
        self.is_wgs = None

        self.min_allele_fraction = None
        self.coverage_interval = None
        self.variant_regions_bed = None
        self.coverage_bed = None
        self.parent_project = parent_project

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


class BaseBatch:
    def __init__(self, name=None, normal=None, tumor=None, rna_samples=None,
                 parent_project=None, somatic_vcf=None, germline_vcf=None, sv_vcf=None):
        self.name = name
        self.parent_project = parent_project

        self.normal = normal
        self.normals = []
        if isinstance(normal, Iterable):
            self.normals = list(normal)
        elif normal:
            self.normals = [normal]

        self.tumor = tumor
        self.tumors = []
        if isinstance(tumor, Iterable):
            self.tumors = list(tumor)
        elif tumor:
            self.tumors = [tumor]

        self.rna_samples = list(rna_samples) if rna_samples else []

        self.somatic_vcf = somatic_vcf
        self.germline_vcf = somatic_vcf
        self.sv_vcf = germline_vcf

        self.somatic_caller = None
        self.germline_caller = None
        self.sv_caller = None

    def is_paired(self):
        return self.normal and self.tumor

    def is_germline(self):
        return self.tumor and self.tumor.phenotype == 'germline'

    def __str__(self):
        return self.name

    def find_somatic_vcf(self, silent=False):
        pass

    def find_germline_vcf(self, silent=False):
        pass

    def find_sv_vcf(self, silent=False):
        pass

    def add_normal(self, sample):
        self.normals.append(sample)
        if not self.normal:
            self.normal = sample

    def add_tumor(self, sample):
        self.tumors.append(sample)
        if not self.tumor:
            self.tumor = sample

    def add_rna_sample(self, sample):
        self.rna_samples.append(sample)


# class Caller:
#     def __init__(self, name=None, is_germline=False):
#         self.name = name
#         self.is_germline = is_germline
#         self.samples = []
#
#     def __str__(self):
#         return self.name