import glob
import re
from collections import defaultdict
import yaml
from os import listdir
from os.path import join, abspath, pardir, splitext, basename, dirname, realpath, isdir, isfile, exists

from ngs_utils.file_utils import adjust_path, verify_dir, file_exists, safe_mkdir, verify_file, add_suffix
from ngs_utils.logger import critical, debug, info, err, warn
from ngs_utils.Sample import BaseSample, BaseBatch, BaseProject
from natsort import natsort_keygen


class DragenProject(BaseProject):
    class DragenSample(BaseSample):
        natsort_key = natsort_keygen()

        def __init__(self, **kwargs):
            BaseSample.__init__(self, **kwargs)  # name, dirpath, work_dir, bam, vcf, phenotype, normal_match
            self.qc_files = []
            self.bam = None

    class DragenBatch(BaseBatch):
        def __init__(self, **kwargs):
            BaseBatch.__init__(self, **kwargs)
            self.batch_qc_files = []
            self.somatic_vcf = join(self.parent_project.dir, self.name + '.hard-filtered.vcf.gz')
            self.sv_vcf = join(self.parent_project.dir, self.name + '.sv.vcf.gz')
            self.replay_file = join(self.parent_project.dir, self.name + '-replay.json')

        def find_somatic_vcf(self, silent=False):
            if isfile(self.somatic_vcf):
                verify_file(self.somatic_vcf, is_critical=True)
                if not silent:
                    info(f'Found somatic VCF in <dragen-dir>/<tumor-name>-hard-filtered.vcf.gz: ' + self.somatic_vcf)

        def find_sv_vcf(self, silent=False):
            if isfile(self.sv_vcf):
                verify_file(self.sv_vcf, is_critical=True)
                if not silent:
                    info(f'Found SV VCF in <dragen-dir>/<tumor-name>.sv.vcf.gz: ' + self.sv_vcf)

        def all_qc_files(self):
            return self.batch_qc_files + self.tumor.qc_files + self.normal.qc_files

        def add_tumor(self, name):
            sample = DragenProject.DragenSample(name=name, phenotype='tumor', batch=self)
            sample.bam = join(self.parent_project.dir, self.name + '_tumor.bam')
            self.tumor = sample
            self.parent_project.samples.append(sample)
            return sample

        def add_normal(self, name):
            sample = DragenProject.DragenSample(name=name, phenotype='normal', batch=self)
            sample.bam = join(self.parent_project.dir, self.name + '.bam')
            self.normal = sample
            self.parent_project.samples.append(sample)
            return sample

    def __init__(self, input_dir=None, silent=False, exclude_samples=None, **kwargs):
        BaseProject.__init__(self, input_dir=input_dir, **kwargs)
        self.somatic_caller = 'dragen'
        self.germline_caller = 'dragen'
        self.genome_build = 'hg38'

        self.bam_list_csv = join(self.dir, 'bam_list.csv')
        debug(f'Parsing project {input_dir}')
        for replay_file in glob.glob(join(self.dir, '*-replay.json')):
            batch_name = basename(replay_file.split('-replay.json')[0])
            debug(f'Found somatic variants for batch {batch_name}')
            if exclude_samples and batch_name in exclude_samples:
                continue
            batch = self.add_batch(batch_name)

            # Reading bam_list.csv to get the sample names. Typical DRAGEN usage includes setting
            # the parameters --RGID and --RGSM-tumor, e.g.:
            # --RGID P025_N --RGSM P025_N --RGID-tumor P025_T --RGSM-tumor P025_T
            # --output-directory /output/P025 --output-file-prefix P025
            # Those values are used inside the VCF files:
            # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  P025_N  P025_T
            # And the VCF files themselves, like all other output files, are prefixed with the --output-file-prefix
            # parameter e.g. P025.hard-filtered.vcf.gz and P025.sv.vcf.gz. The only way to get the RGSM values
            # is from the bam_list.csv file which is not prefixed with P025 and thus can be accidentally
            # overriden. Taking this risk and reading that file:
            # RGID,SampleID,Library,Lane,BamFile
            # P025_N,0x55d4760,0x55d87b0,0x55d87d0,/output/P025/P025.bam
            # P025_T,0x55d4768,0x55d87b8,0x55d87d8,/output/P025/P025_tumor.bam
            tumor_name = None
            normal_name = None
            with open(self.bam_list_csv) as f:
                for i, l in enumerate(f):
                    if i > 0:
                        sn, _, _, _, bam_path = l.strip().split(',')
                        if basename(bam_path) == batch_name + '_tumor.bam':
                            tumor_name = sn
                        if basename(bam_path) == batch_name + '.bam':
                            normal_name = sn
            assert tumor_name, f'Cannot find tumor sample name in {self.bam_list_csv}'
            assert normal_name, f'Cannot find normal sample name in {self.bam_list_csv}'

            batch.add_tumor(tumor_name)
            batch.add_normal(normal_name)
            if exclude_samples and batch.normal.name in exclude_samples:
                continue
            self.samples.extend([batch.tumor, batch.normal])
            batch.tumor.bam = join(self.dir, batch_name + '_tumor.bam')
            batch.normal.bam = join(self.dir, batch_name + '.bam')
            self.batch_by_name[batch_name] = batch

            batch.find_somatic_vcf(silent=silent)
            batch.find_germline_vcf(silent=silent)
            batch.find_sv_vcf(silent=silent)

            # populating qc files for multiqc:
            for suffix in [
                '.fragment_length_hist.csv',
                '.mapping_metrics.csv',
                '.ploidy_estimation_metrics.csv',
                '.time_metrics.csv',
                '.vc_metrics.csv',
                '.sv_metrics.csv',
            ]:
                qc_fpath = join(self.dir, f'{batch_name}{suffix}')
                if isfile(qc_fpath):
                    debug(f'Found QC file for batch {qc_fpath}')
                    batch.batch_qc_files.append(qc_fpath)
                else:
                    debug(f'Can\'t find QC file for batch {qc_fpath}')

            for suffix in [
                '.wgs_contig_mean_cov_{phenotype}.csv',
                '.wgs_coverage_metrics_{phenotype}.csv',
                '.wgs_fine_hist_{phenotype}.csv',
            ]:
                qc_fpath = join(self.dir, f'{batch_name}{suffix}')
                if isfile(qc_fpath.format(phenotype="normal")):
                    debug(f'Found QC file for normal sample {qc_fpath}')
                    batch.normal.qc_files.append(qc_fpath.format(phenotype="normal"))
                else:
                    debug(f'Can\'t find QC file for normal sample {qc_fpath}')
                if isfile(qc_fpath.format(phenotype="tumor")):
                    debug(f'Found QC file for tumor sample {qc_fpath}')
                    batch.tumor.qc_files.append(qc_fpath.format(phenotype="tumor"))
                else:
                    debug(f'Can\'t find QC file for tumor sample {qc_fpath}')

            debug(f'Found {len(batch.batch_qc_files)} batch QC files, '
                  f'{len(batch.tumor.qc_files)} tumor QC files, '
                  f'{len(batch.normal.qc_files)} normal QC files')

        if len(self.batch_by_name) == 1:
            self.project_name = list(self.batch_by_name.values())[0].name
        else:
            self.project_name = basename(input_dir)

    def add_batch(self, batch_name):
        batch = DragenProject.DragenBatch(name=batch_name, parent_project=self)
        self.batch_by_name[batch_name] = batch
        return batch

