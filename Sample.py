class BaseSample:
    def __init__(self, name, dirpath, bam=None, bed=None, vcf=None, genome=None,
                 targqc_dirpath=None, fastqc_dirpath=None, picard_dirpath=None, clinical_report_dirpath=None,
                 flagged_regions_dirpath=None, normal_match=None, sv_fpath=None, sv_bed=None):
        self.name = name
        self.bam = bam
        self.dedup_bam = None
        self.bed = bed
        self.sv_bed = sv_bed
        self.is_wgs = False
        self.vcf = vcf
        self.dirpath = dirpath
        self.phenotype = None
        self.gender = None
        self.genome = None
        self.var_dirpath = None
        self.normal_match = normal_match
        self.min_af = None
        self.sv_fpath = sv_fpath
        self.targqc_dirpath = targqc_dirpath

    def __cmp__(self, other):
        return cmp(self.key_to_sort(), other.key_to_sort())

    def key_to_sort(self):
        parts = []

        cur_part = []
        prev_was_num = False

        for c in self.name:
            if prev_was_num == c.isdigit() and c not in ['-', '.']:  # same type of symbol, but not - or .
                cur_part.append(c)
            else:
                if cur_part:
                    part = ''.join(cur_part)
                    if prev_was_num:
                        part = int(part)
                    parts.append(part)
                    cur_part = []

                if c in ['-', '.']:
                    pass
                else:
                    if c.isdigit():
                        prev_was_num = True
                    else:
                        prev_was_num = False
                    cur_part.append(c)
        if cur_part:
            part = ''.join(cur_part)
            if prev_was_num:
                part = int(part)
            parts.append(part)

        return tuple(parts)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    @classmethod
    def load(cls, data):
        sample = cls(**data)
        sample.__dict__ = data
        return sample
