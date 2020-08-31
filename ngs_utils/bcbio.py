import re
import shutil
import tarfile
from collections import defaultdict
import yaml
from os import listdir
from os.path import join, abspath, pardir, splitext, basename, dirname, realpath, isdir, isfile, exists

from ngs_utils.Sample import BaseSample, BaseBatch, BaseProject
from ngs_utils.bam_utils import verify_bam
from ngs_utils.call_process import run, run_simple
from ngs_utils.config import load_yaml_config
from ngs_utils.file_utils import adjust_path, verify_dir, file_exists, safe_mkdir, verify_file, add_suffix
from ngs_utils.logger import critical, debug, info, err, warn
from ngs_utils.key_genes_utils import get_target_genes, is_small_target


CALLER_PRIORITY = ['ensemble', 'strelka2', 'vardict', 'gatk-haplotype']


class BcbioSample(BaseSample):
    def __init__(self, **kwargs):
        BaseSample.__init__(self, **kwargs)
        self.old_name = None
        self.raw_name = None
        self.project_tag = None
        self.genome_build = None
        self.sample_info = dict()
        self.sv_regions_bed = None
        self.variantcallers = []

        self.germline_caller = None
        self.somatic_caller = None

    def get_name_for_files(self):  # In case if the sample if symlink from another project, and the name was changed in this one
        return self.old_name or self.name

    @staticmethod
    def parse_sample_ids(sample_info):
        description = str(sample_info['description']).replace('.', '_')

        batch_names = sample_info.get('metadata', dict()).get('batch')
        if isinstance(batch_names, int) or isinstance(batch_names, float):
            batch_names = str(batch_names)
        if isinstance(batch_names, str):
            batch_names = [batch_names]
        if batch_names:
            batch_names = [b.replace('.', '_') for b in batch_names if b]

        return description, batch_names

    @staticmethod
    def load_from_sample_info(sample_info, bcbio_project,
                              include_samples=None, exclude_samples=None,
                              extra_batches=None, silent=False):
        """ Get sample and batch names and exclude/include based on exclude_samples and include_samples
        """
        description, batch_names = BcbioSample.parse_sample_ids(sample_info)

        if exclude_samples:
            # Sample name
            if description in exclude_samples:
                if not silent: info(f'Skipping sample {description}')
                return None
            # Batch names
            if batch_names:
                filtered_batch_names = [b for b in batch_names if b not in exclude_samples]
                if not filtered_batch_names:
                    if not silent: info(f'Skipping sample {description} with batch info {", ".join(batch_names)}')
                    return None
                batch_names = filtered_batch_names

        if include_samples:
            # Sample name
            if description in include_samples:
                if not silent: info(f'Using sample {description} and all samples sharing batches {batch_names}')
            else:
                # Batch names
                if batch_names:
                    incl_batch_names = [b for b in batch_names if b in include_samples]
                    if incl_batch_names:
                        if not silent: info(f'Using sample {description} with batch info {", ".join(batch_names)}')
                    extr_batch_names = [b for b in batch_names if extra_batches and b in extra_batches]
                    if extr_batch_names and not incl_batch_names:
                        if not silent: info(f'Using sample {description} as it shares batches {extr_batch_names} with included samples')
                    incl_batch_names += extr_batch_names

                    if incl_batch_names:
                        batch_names = incl_batch_names
                    else:
                        return None

        # Creating BcbioSample object
        s = BcbioSample(parent_project=bcbio_project)
        s.sample_info = sample_info
        if 'description_original' in sample_info:
            s.old_name = str(sample_info['description_original']).replace('.', '_')

        # Setting phenotype and batches
        s.phenotype = sample_info.get('metadata', dict()).get('phenotype', 'tumor')
        if not batch_names:
            batch_names = [s.get_name_for_files() + '-batch']
        if len(batch_names) > 1 and s.phenotype != 'normal':
            critical('Multiple batches for non-normal ' + s.phenotype + ' sample ' + s.name + ': ' + ', '.join(batch_names))
        s.batch_names = batch_names

        # Setting genome build based reference paths
        s.genome_build = sample_info['genome_build']
        s.variant_regions_bed = s.parent_project.config_path(val=sample_info['algorithm'].get('variant_regions'))
        s.sv_regions_bed = s.parent_project.config_path(val=sample_info['algorithm'].get('sv_regions')) or s.variant_regions_bed
        s.coverage_bed = s.parent_project.config_path(val=sample_info['algorithm'].get('coverage')) or s.sv_regions_bed
        if s.coverage_bed and not isfile(s.coverage_bed):
            if not silent:
                debug('coverage bed ' + str(s.coverage_bed) + ' not found. Looking relatively to genomes "basedir"')
            try:
                import az
            except ImportError:
                pass
            else:
                genome_cfg = az.get_refdata(s.genome_build)
                ref_basedir = genome_cfg.get('basedir')
                if not ref_basedir:
                    critical('coverage bed ' + str(s.coverage_bed) + ' not found and "basedir" not provided in system config')
                s.coverage_bed = join(ref_basedir, 'coverage', 'prioritize', s.coverage_bed) + '.bed'

        s.is_rnaseq = 'rna' in sample_info['analysis'].lower()
        s.min_allele_fraction = (1.0/100) * float(sample_info['algorithm'].get('min_allele_fraction', 1.0))
        if s.variant_regions_bed is None:
            s.coverage_interval = 'genome'
        else:
            s.coverage_interval = 'regional'
        s.is_wgs = s.coverage_interval == 'genome'

        s._set_name_and_paths(
            name=description,
            variantcallers_data=sample_info['algorithm'].get('variantcaller'),
            ensemble='ensemble' in sample_info['algorithm'],
            silent=silent)
        return s

    def find_bam(self, silent=False):
        name = self.get_name_for_files()

        to_try = [
            '-ready.bam',
            '-ready.cram',
            '-sort.bam',
        ]
        for ext in to_try:
            fpath = adjust_path(join(self.dirpath, name + ext))
            if verify_file(fpath):
                return fpath

        input_file = self.sample_info['files']
        if not isinstance(input_file, str):
            input_file = input_file[0]
        if isinstance(input_file, str) and input_file.endswith('.bam'):
            debug('Bcbio was run from BAM input')
            if not input_file.startswith('/'):
                input_file = abspath(join(self.parent_project.work_dir, input_file))
            if verify_file(input_file):
                debug('Using BAM file from input YAML ' + input_file)
                return input_file
            else:
                debug('Input BAM file for sample ' + self.name + ' in YAML ' + input_file + ' does not exist')

        if not silent:
            warn('No BAM or CRAM file found for ' + self.name)

    def _set_name_and_paths(self, name, variantcallers_data, ensemble=False, silent=False):
        self.raw_name = name
        self.name = self.raw_name.replace('.', '_')
        self.rgid = self.name
        self.dirpath = verify_dir(join(self.parent_project.final_dir, self.name))
        if not verify_dir(self.dirpath, silent=silent):
            critical(f'Sample "{self.name}" specified in bcbio YAML is not found in the final directory '
                     f'{self.parent_project.final_dir}. Please check consistency between the YAML '
                     f'{self.parent_project.bcbio_yaml_fpath} and the directories in `final`: '
                     f'to every "description" value in YAML, there should be a corresponding folder with the '
                     f'same name in `final`. You can use `-e` option to exclude samples (comma-separated) '
                     f'from consideration, if you are sure that missing folders are expected.')

        self.bam = self.find_bam(silent=silent)

        if self.is_rnaseq:
            gene_counts = adjust_path(join(self.dirpath, self.get_name_for_files() + '-ready.counts'))
            if isfile(gene_counts) and verify_file(gene_counts):
                self.counts_file = gene_counts
            else:
                if not silent: warn('Counts for ' + self.name + ' not found')
        else:
            if variantcallers_data:
                self._set_variant_callers(variantcallers_data, ensemble=ensemble)
            else:
                if not silent: warn('No variant callers set in config, skipping finding VCF files')

    def _set_variant_callers(self, variantcallers_data, ensemble=False):
        if isinstance(variantcallers_data, dict):
            if 'germline' in variantcallers_data and self.phenotype == 'normal':
                self.variantcallers = variantcallers_data.get('germline')
            else:
                self.variantcallers = variantcallers_data.get('somatic')

        if isinstance(variantcallers_data, str):
            self.variantcallers = [variantcallers_data]
        elif isinstance(variantcallers_data, list):
            self.variantcallers = variantcallers_data

        if ensemble and len(self.variantcallers) > 1:
            self.variantcallers = ['ensemble'] + self.variantcallers

        if self.phenotype != 'germline' and self.phenotype != 'normal':
            self.somatic_caller = next((c for c in CALLER_PRIORITY if c in self.variantcallers),
                                       self.variantcallers[0])
        else:
            self.germline_caller = next((c for c in CALLER_PRIORITY if c in self.variantcallers),
                                        self.variantcallers[0])

    def find_sv_vcf(self):
        return self.find_cnv_file(self.name + '-manta.vcf.gz') or \
               self.find_cnv_file(self.name + '-lumpy.vcf.gz')

    def find_cnv_file(self, fname):
        for fpath in [join(self.dirpath, fname)]:
            if isfile(fpath):
                return verify_file(fpath, silent=True)

    def find_coverage_stats(self):
        sname = self.name
        dirpath = self.dirpath
        if self.phenotype == 'germline':
            sname = re.sub(r'-germline$', '', sname)
            dirpath = re.sub(r'-germline$', '', dirpath)
        return verify_file(join(dirpath, 'qc', 'coverage', sname + '_coverage.bed'), silent=True)

    def get_metric(self, names):
        if isinstance(names, str):
            names = [names]
        if not self.sample_info or not self.sample_info.get('metrics'):
            return None
        metrics = self.sample_info['metrics']
        val = None
        for k in metrics:
            if k.lower() in [n.lower() for n in names] and metrics[k] != 'NA':
                val = metrics[k]
        if val is None:
            err('Cannot find ' + ', '.join(names) + ' in metrics for ' + self.name)
        return val

    def get_avg_depth(self):
        return self.get_metric(['Avg_coverage', 'Avg_coverage_per_region'])

    def get_reads_count(self):
        return self.get_metric(['Total_reads', 'Total reads'])

    def get_usable_count(self):
        if self.get_metric('Usable_pct'):
            return int(self.get_reads_count() * self.get_metric('Usable_pct') / 100)

    def is_dedupped(self):
        return self.sample_info.get('algorithm', {}).get('mark_duplicates', False)


class BcbioBatch(BaseBatch):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.somatic_caller = 'ensemble'
        self.germline_caller = 'ensemble'
        self.sv_caller = 'manta'

    def find_somatic_vcf(self, silent=False, caller=None):
        caller = caller or self.somatic_caller
        if not caller:
            if not silent:
                warn(f'Batch {self.name} have no variant caler info assigned, skipping finding somatic VCF')
                return

        # in datestamp. cwl-bcbio writes there
        vcf_cwl_fpath_gz = adjust_path(join(self.parent_project.date_dir, self.name + '-' + caller + '.vcf.gz'))
        # in datestamp. bcbio before 1.1.6
        vcf_old_fpath_gz = adjust_path(join(self.parent_project.date_dir, self.name + '-' + caller + '-annotated.vcf.gz'))
        # in sample dir. starting from bcbio 1.1.6, ~ Dec 2019
        vcf_fpath_gz = adjust_path(join(self.tumor.dirpath, self.tumor.name + '-' + caller + '.vcf.gz'))

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
        caller = caller or self.germline_caller
        if not caller:
            if not silent:
                warn(f'Batch {self.name} have no variant caler info assigned, skipping finding germline VCF')
            return
        assert caller

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

    def find_sv_vcf(self, silent=False, caller=False):
        caller = caller or self.sv_caller

        sv_prio   = join(self.tumor.dirpath, f'{self.name}-sv-prioritize-{caller}.vcf.gz')
        sv_unprio = join(self.tumor.dirpath, f'{self.name}-{caller}.vcf.gz')
        # CWL?
        sv_cwl_prio   = join(self.parent_project.date_dir, f'{self.tumor.name}-{caller}-prioritized.vcf.gz')
        sv_cwl_unprio = join(self.parent_project.date_dir, f'{self.tumor.name}-{caller}.vcf.gz')

        if isfile(sv_prio):
            verify_file(sv_prio, is_critical=True)
            if not silent: info(f'Found SV VCF in <tumor>/<batch>-sv-prioritize-{caller}.vcf.gz: ' + sv_prio)
            self.sv_vcf = sv_prio

        elif isfile(sv_unprio):
            verify_file(sv_unprio, is_critical=True)
            if not silent: info(f'Found SV VCF in <tumor>/<batch>-{caller}.vcf.gz: ' + sv_unprio)
            self.sv_vcf = sv_unprio

        elif isfile(sv_cwl_prio):
            verify_file(sv_cwl_prio, is_critical=True)
            if not silent: info(f'Found SV VCF in <date-dir>/<tumor-name>-{caller}-prioritized.vcf.gz: ' + sv_cwl_prio)
            self.sv_cwl_prio = sv_cwl_prio

        elif isfile(sv_cwl_unprio):
            verify_file(sv_cwl_unprio, is_critical=True)
            if not silent: info(f'Found SV VCF in <date-dir>/<tumor-name>-{caller}.vcf.gz: ' + sv_cwl_prio)
            self.sv_vcf = sv_cwl_unprio

        elif not silent:
            warn(f'Could not find SV VCF file for batch {self.name}, caller {caller} neither under sample folder as '
                 f'<tumor>/<batch>(-sv-prioritize)-{caller}.vcf.gz (conventional bcbio), '
                 f'nor in the project folder as project/<tumor>-{caller}(-prioritized).vcf.gz (CWL bcbio).')

    def find_qc_files(self, dst_dir, exclude_files=None, include_files=None):
        """
        Parses bcbio MultiQC file list and collects all QC files belonging to this batch

        :param dst_dir: destination directory where the QC files will be copied to
        :param exclude_files: not include files matching these patterns
        :param include_files: only include files matching these patterns
        :return: list of file paths copied into `new_mq_data_dir`
        """

        mq_dir = join(self.parent_project.date_dir, 'multiqc')
        mq_filelist = join(mq_dir, 'list_files_final.txt')
        verify_file(mq_filelist, is_critical=True)

        # Cromwell?
        cwl_targz = join(mq_dir, 'multiqc-inputs.tar.gz')
        tar_f_by_fp = dict()
        if isfile(cwl_targz):
            info(f'Found CWL MultiQC output {cwl_targz}, extracting required QC files from the archive')
            if cwl_targz:
                tar = tarfile.open(cwl_targz)
                for member in tar.getmembers():
                    rel_fp = member.name
                    if 'call-multiqc_summary/execution/qc/multiqc/' in rel_fp:
                        rel_fp = rel_fp.split('call-multiqc_summary/execution/qc/multiqc/')[1]
                    tar_f_by_fp[rel_fp] = tar.extractfile(member)

        qc_files_not_found = []
        qc_files_found = []
        with open(mq_filelist) as inp:
            for fp in [l.strip() for l in inp if l.strip()]:
                if fp == 'trimmed' or fp.endswith('/trimmed'):
                    continue  # back-compatibility with bcbio
                if exclude_files:
                    if isinstance(exclude_files, str):
                        exclude_files = [exclude_files]
                    if any(re.search(ptn, fp) for ptn in exclude_files):
                        continue
                if include_files:
                    if isinstance(include_files, str):
                        include_files = [include_files]
                    if not any(re.search(ptn, fp) for ptn in include_files):
                        continue

                new_fp = _extract_qc_file(fp, dst_dir, self.parent_project.final_dir, tar_f_by_fp)
                if not new_fp:
                    qc_files_not_found.append(fp)
                    continue
                else:
                    qc_files_found.append(new_fp)

        if qc_files_not_found:
            warn('-')
            warn(f'Some QC files from list {mq_filelist} were not found:' +
                ''.join('\n  ' + fpath for fpath in qc_files_not_found))
        return qc_files_found


class NoConfigDirException(Exception):
    pass
class NoDateStampsException(Exception):
    pass
class MultipleDateStampsException(Exception):
    pass


class BcbioProject(BaseProject):
    def __init__(self, input_dir=None, project_name=None, proc_name='postproc',
                 include_samples=None, exclude_samples=None, silent=False, **kwargs):
        super().__init__(input_dir, **kwargs)
        self.config_dir = None
        self.final_dir = None
        self.date_dir = None
        self.log_dir = None
        self.postproc_log_dir = None
        self.work_dir = None
        self.bcbio_yaml_fpath = None

        self.versions = None
        self.programs = None

        self.samples_by_caller = defaultdict(list)  # (caller, is_germline) -> [samples]

        self.variant_regions_bed = None
        self.sv_regions_bed = None          # "sv_regions" or "variant_regions"
        self.original_coverage_bed = None
        self.coverage_bed = None            # "coverage" or "sv_regions" or "variant_regions"

        self.coverage_interval = None  # amplicon, regional, genome
        self.min_allele_fraction = None
        self.is_wgs = None
        self.is_rnaseq = None
        self.postproc_mqc_files = []

        self.silent = silent

        if input_dir:
            self.load_from_bcbio_dir(
                input_dir, project_name, proc_name,
                include_samples=include_samples, exclude_samples=exclude_samples)

    def set_project_level_dirs(self, bcbio_cnf, config_dir, project_name=None, final_dir=None, date_dir=None,
                               create_dirs=False, proc_name='postproc'):
        self.final_dir = self.set_final_dir(bcbio_cnf, config_dir, final_dir)
        if create_dirs: safe_mkdir(self.final_dir)

        self.project_name = self._set_project_name(self.final_dir, project_name)

        self.work_dir = abspath(join(self.final_dir, pardir, 'work'))
        if create_dirs: safe_mkdir(self.work_dir)

        self.date_dir = self._set_date_dir(bcbio_cnf, self.final_dir, date_dir, create_dir=create_dirs,
                                           silent=self.silent)
        self.log_dir = join(self.date_dir, 'log')
        self.postproc_log_dir = join(self.log_dir, proc_name)
        if create_dirs: safe_mkdir(self.postproc_log_dir)

        self.versions = verify_file(join(self.date_dir, 'data_versions.txt'), silent=True)
        self.programs = verify_file(join(self.date_dir, 'programs.txt'), silent=True)

    def load_from_bcbio_dir(self, input_dir, project_name=None, proc_name='postproc',
                            include_samples=None, exclude_samples=None):
        """
        Analyses existing bcbio folder.
        input_dir: root bcbio folder, or any other directory inside it
        """
        self.config_dir, self.final_dir, self.date_dir = detect_bcbio_dir(input_dir, silent=self.silent)
        bcbio_cnf, self.bcbio_yaml_fpath = load_bcbio_cnf(self.config_dir, silent=self.silent)
        self.set_project_level_dirs(bcbio_cnf, self.config_dir, project_name=project_name, final_dir=self.final_dir,
                                    date_dir=self.date_dir, proc_name=proc_name)
        self.set_samples(bcbio_cnf, include_samples=include_samples, exclude_samples=exclude_samples)
        self._load_bcbio_summary()
        # self._load_target_info()
        return self

    def set_samples(self, bcbio_cnf, include_samples=None, exclude_samples=None):
        debug('Reading sample details...')
        exclude_samples = [s.replace('.', '_') for s in exclude_samples] if exclude_samples else None
        include_samples = [s.replace('.', '_') for s in include_samples] if include_samples else None

        # First pass - just to get extra batch IDs that we need to include to have batches consistent
        extra_batches = set()
        all_sample_names = set()
        all_batch_names = set()
        if include_samples:
            for sample_info in bcbio_cnf['details']:
                sname, batch_names = BcbioSample.parse_sample_ids(sample_info)
                all_sample_names.add(sname)
                all_batch_names |= set(batch_names)
                if sname in include_samples:
                    for b in batch_names:
                        if b not in (include_samples or []) and b not in (exclude_samples or []):
                            extra_batches.add(b)

        # Second pass - including/excluding, and creating BcbioSample objects
        for sample_info in bcbio_cnf['details']:
            s = BcbioSample.load_from_sample_info(
                sample_info,
                bcbio_project=self,
                include_samples=include_samples,
                exclude_samples=exclude_samples,
                extra_batches=extra_batches,
                silent=self.silent)
            if s:
                self.samples.append(s)

        if not self.samples:
            if exclude_samples:
                critical(f'Error: no samples left with the exclusion of '
                         f'batch/sample name(s): {", ".join(exclude_samples)}\n'
                         f'Available samples from the YAML file {self.bcbio_yaml_fpath}:\n'
                         f'{", ".join(all_sample_names)}\nbatches: {", ".join(all_batch_names)}')
            if include_samples:
                critical(f'Error: could not find a batch or a sample with the name(s): '
                         f'{", ".join(include_samples)}\n'
                         f'Available samples from the YAML file {self.bcbio_yaml_fpath}:\n'
                         f'{", ".join(all_sample_names)}\nbatches: {", ".join(all_batch_names)}')
            critical(f'Error: could not parse any batch or samples in the bcbio project. '
                     f'Please check the bcbio YAML file: {self.bcbio_yaml_fpath}')

        not_found_samples = [s.name for s in self.samples if not s.bam]
        if not_found_samples:
            if not self.silent: warn(f'Warning: no BAM files not found for {len(not_found_samples)}/{len(self.samples)} samples')

        self.samples.sort(key=lambda _s: _s.key_to_sort())
        self.batch_by_name = self.update_batches(self.samples, self.silent)

        def _check_dup_props(prop, is_critical=False):
            _vals = set([s_.__dict__.get(prop) for s_ in self.samples])
            if len(_vals) > 1:
                (critical if is_critical else err)('Got different ' + prop + ' values in samples in ' + self.project_name)
            else:
                self.__dict__[prop] = _vals.pop()
        _check_dup_props('genome_build')
        _check_dup_props('variant_regions_bed')
        _check_dup_props('coverage_bed')
        _check_dup_props('sv_regions_bed')
        _check_dup_props('is_rnaseq')
        _check_dup_props('min_allele_fraction')
        _check_dup_props('is_wgs', is_critical=False)
        _check_dup_props('coverage_interval', is_critical=False)
        if self.is_rnaseq:
            debug('RNAseq')
        elif self.coverage_interval:
            debug('Coverage interval: ' + str(self.coverage_interval))

        for s in self.samples:
            for caller in s.variantcallers:
                self.samples_by_caller[(caller, s.phenotype == 'germline')].append(s)

        debug('Done loading bcbio project ' + self.project_name)

    def _load_bcbio_summary(self):
        fp = self.find_in_log('project-summary.yaml')
        if fp:
            with open(fp) as f:
                data = yaml.load(f)
            metrics_by_sample = dict()
            for s_data in data.get('samples', []):
                metrics_by_sample[s_data['description']] = s_data.get('summary', dict()).get('metrics')
            for s in self.samples:
                sname = s.name
                if s.phenotype == 'germline':
                    sname = re.sub(r'-germline$', '', s.name)
                s.sample_info['metrics'] = metrics_by_sample[sname]

    def config_path(self, val):
        if not val:
            return val
        full_path = adjust_path(join(self.config_dir, val))
        if exists(full_path):
            return full_path
        else:
            return val

    @staticmethod
    def _set_date_dir(bcbio_cnf, final_dir, date_dir, create_dir=False, silent=False):
        if not date_dir:
            fc_date = bcbio_cnf.get('fc_date')
            fc_name = bcbio_cnf.get('fc_name') or 'project'
            if fc_date:
                # Date dirpath is from bcbio and named after fc_name, not our own project name
                date_dir = join(final_dir, fc_date + '_' + fc_name)
                if not create_dir and not verify_dir(date_dir, silent=True):
                    critical('Error: no project directory of format {fc_date}_{fc_name} or {fc_name}_{fc_date}')
            else:
                if isdir(join(final_dir, 'project')):  # bcbio-CWL?
                    date_dir = join(final_dir, 'project')
                    if not silent: info('Using the datestamp dir from bcbio-CWL: ' + date_dir)
                else:
                    regexs = [fr'^\d\d\d\d-[01][0-9]-[0-3][0-9]_{fc_name}']
                    date_dirs = [join(final_dir, dirpath)
                                 for dirpath in listdir(final_dir)
                                 if any(re.match(regex, dirpath) for regex in regexs)]
                    if len(date_dirs) == 0:
                        raise NoDateStampsException('Error: no datestamp directory!')
                    elif len(date_dirs) == 1:
                        date_dir = date_dirs[0]
                    else:
                        dates = [(tuple(map(int, basename(d).split('_')[0].split('-'))), d) for d in date_dirs]
                        newest_date, newest_dir = sorted(dates, reverse=True)[0]
                        newest_dirs = [d_dir for d_dir in date_dirs if d_dir == newest_dir]
                        if len(newest_dirs) > 1:
                            raise MultipleDateStampsException(f'Error: multiple datestamp directory found, '
                               f'and can\'t select the most recent one because there are multiple latest dirs: {newest_dirs}')
                        date_dir = newest_dirs[0]

                    if not silent: info('Using the datestamp dir: ' + date_dir)
        if create_dir:
            safe_mkdir(date_dir)
        return date_dir

    @staticmethod
    def _set_project_name(final_dir, project_name=None):
        if not project_name:
            root_dir = dirname(final_dir)
            # path is like ../Bio_0031_Heme_MRL_DLBCL_IRAK4/bcbio_Dev_0079/final
            second_part = basename(root_dir)  # bcbio_Dev_0079
            bcbio_project_parent_dirname = basename(dirname(root_dir))  # Bio_0031_Heme_MRL_DLBCL_IRAK4
            project_name = bcbio_project_parent_dirname + '_' + second_part
        return project_name

    @staticmethod
    def set_final_dir(bcbio_cnf, config_dir, final_dir=None, create_dir=False):
        if final_dir:
            return final_dir
        elif 'upload' in bcbio_cnf and 'dir' in bcbio_cnf['upload']:
            final_dirname = bcbio_cnf['upload']['dir']
            final_dir = adjust_path(join(config_dir, final_dirname))
            if create_dir: safe_mkdir(final_dir)
            verify_dir(final_dir, 'upload directory specified in the bcbio config', is_critical=True)
        else:
            final_dir = abspath(join(config_dir, pardir, 'final'))
            if create_dir: safe_mkdir(final_dir)
            if not verify_dir(final_dir):
                critical('If final directory it is not named "final", please, specify it in the bcbio config.')
        return final_dir

    def update_batches(self, samples, silent=False):
        batch_by_name = {bn: BcbioBatch(name=bn, parent_project=self)
                         for bn in list(set([b for s in samples for b in s.batch_names]))}

        for sample in samples:
            for bn in sample.batch_names:
                batch_by_name[bn].name = bn
                sample.batches.append(batch_by_name[bn])
                if sample.phenotype == 'normal':
                    if batch_by_name[bn].normal:
                        critical('Multiple normal samples for batch ' + bn)
                    batch_by_name[bn].add_normal(sample)
                else:
                    batch_by_name[bn].add_tumor(sample)

        # Removing batches that do not have matching tumor samples
        batch_by_name = {bn: b for bn, b in batch_by_name.items() if b.tumor}

        # for batch in batch_by_name.values():
        #     if batch.normal and not batch.tumor:
        #         if not silent: info('Batch ' + batch.name + ' contains only normal, treating sample ' + batch.normal.name + ' as tumor')
        #         batch.normal.phenotype = 'tumor'
        #         batch.normal.batch = batch
        #         batch.tumor = batch.normal
        #         batch.normal = None

        # setting up batch properties
        for b in batch_by_name.values():
            if b.tumor:
                b.tumor.normal_match = b.normal

        # setting variant caller names for batches
        for b in batch_by_name.values():
            if b.tumor.somatic_caller is None:
                if not silent:
                    warn(f'Sample {b.tumor} doesn\'t have somatic variant callers info, skip assinging '
                         f'variant caller to batch {b.name}')
            else:
                b.somatic_caller = b.tumor.somatic_caller
            if b.normal:
                if b.normal.germline_caller is None:
                    if not silent:
                        warn(f'Sample {b.tumor} doesn\'t have germline variant callers info, skip assinging '
                             f'germline variant caller to batch {b.name}')
                else:
                    b.germline_caller = b.normal.germline_caller

        # finding vcfs
        for b in batch_by_name.values():
            if b.tumor:
                b.find_somatic_vcf(silent=silent)
                b.find_sv_vcf(silent=silent)
            if b.normal:
                b.find_germline_vcf(silent=silent)

        return batch_by_name

    def find_in_log(self, fname, is_critical=False, silent=True):
        options = [join(self.log_dir, fname),
                   join(self.date_dir, fname)]
        for fpath in options:
            if isfile(fpath):
                return fpath
        if is_critical:
            critical('Log file not found as ' + ', '.join(options))
        elif not silent:
            err('Log file not found as ' + ', '.join(options))

    def get_target_genes(self, get_key_genes_file=None):
        return get_target_genes(self.genome_build, self.coverage_bed,
                                get_key_genes_file=get_key_genes_file)

    def is_small_target(self):
        return is_small_target(self.coverage_bed)


def _extract_qc_file(fp, new_mq_data_dir, final_dir, f_by_fp=None):
    """ Extracts QC file `fp` either by copying from `final_dir` (native bcbio),
        or from tar.gz file `tar_path` (CWL bcbio). Writes into a new file at new_mq_data_dir
    """
    if fp.startswith('report/metrics/'):
        fp = fp.replace('report/metrics/', 'project/multiqc/')  # for CWL _bcbio.txt files

    dst_fp = join(new_mq_data_dir, fp)

    fp_in_final = join(final_dir, fp)
    if isfile(fp_in_final):
        safe_mkdir(dirname(dst_fp))
        shutil.copy2(fp_in_final, dst_fp)
        return dst_fp

    elif f_by_fp and fp in f_by_fp:
        safe_mkdir(dirname(dst_fp))
        with open(dst_fp, 'wb') as out:
            out.write(f_by_fp[fp].read())
        return dst_fp


def detect_bcbio_dir(input_dir, silent=False):
    """
    :param input_dir: `config` dir, or `final` dir, or datestamp dir, or the directory root to `final`
    :return: (config_dir, final_dir, date_dir)
    """
    config_dir, final_dir, date_dir = None, None, None

    input_dir = abspath(input_dir)

    # We are inside `*final*`
    if 'final' in basename(input_dir):  # allow prefixes and postfixes
        final_dir = input_dir
        root_dir = dirname(final_dir)
        config_dir = join(root_dir, 'config')
        if not isdir(config_dir):
            err(f'Are you running on a bcbio output?\n'
                f'The input folder appear to be `final` ({input_dir}), '
                f'however can\'t find `config` directory at the same level ({config_dir})')
            raise NoConfigDirException('No config dir')

    # We are inside `config`
    elif basename(input_dir) == 'config':
        config_dir = input_dir

    # We are in a parent dir to `config` (and possibly `final`, called otherwise)
    elif isdir(join(input_dir, 'config')):
        config_dir = join(input_dir, 'config')

    # We are inside a date dir
    elif isdir(abspath(join(input_dir, pardir, pardir, 'config'))):
        final_dir = abspath(join(input_dir, pardir))
        root_dir = abspath(join(input_dir, pardir, pardir))
        config_dir = abspath(join(root_dir, 'config'))

        # if 'final' not in basename(final_dir):
        #     err(f'Are you running on a bcbio output?\n'
        #         f'Found config directory 2 level up at {config_dir}, assuming your input {input_dir} '
        #         f'is a datestamp directory. However, the parent directory is not called `*final*`')
        #     raise NoConfigDirException('No final dir')

    else:
        if not silent:
            err(f'Are you running on a bcbio output?\n'
                f'{input_dir} is not `config` or `*final*`, and '
                f'can\'t find a `config` directory at {join(input_dir, "config")}, or {abspath(join(input_dir, pardir, "config"))}.'
                f'Make sure that you changed to a bcbio root or final directory, or provided it as a first argument.')
        raise NoConfigDirException('No config dir')

    if not silent:
        if not silent:
            info(f'Bcbio config directory: ' + config_dir)
        if final_dir:
            if not silent: info('"final" directory: ' + final_dir)
            if date_dir:
                if not silent: info('"datestamp" directory: ' + date_dir)

    return config_dir, final_dir, date_dir


def load_bcbio_cnf(config_dir, silent=False):
    all_yamls = [
        abspath(join(config_dir, fname))
        for fname in listdir(config_dir)
        if fname.endswith('.yaml')]
    if len(all_yamls) == 0:
        critical('No YAML file in the config directory.')

    bcbio_yamls = []
    for fpath in all_yamls:
        if not fpath.endswith('-template.yaml'):
            if 'details' in load_yaml_config(fpath):
                bcbio_yamls.append(fpath)
    if len(bcbio_yamls) == 0:
        critical('No bcbio YAMLs found in the config directory: ' + config_dir +
                 ' (only ' + ', '.join(map(basename, all_yamls)) +
                 ' which do not have the "details" section)')
    if len(bcbio_yamls) > 1:
        critical('More than one bcbio YAML file found in the config directory ' +
                 config_dir + ': ' + ' '.join(bcbio_yamls))
    yaml_fpath = bcbio_yamls[0]
    if not silent: info('Using bcbio YAML config: ' + yaml_fpath)
    return load_yaml_config(yaml_fpath), yaml_fpath


def _normalize(name):
    return name.lower().replace('_', '').replace('-', '')


def ungzip_if_needed(cnf, fpath, silent=False):
    if fpath.endswith('.gz'):
        fpath = fpath[:-3]
    if not file_exists(fpath) and file_exists(fpath + '.gz'):
        gz_fpath = fpath + '.gz'
        cmdline = 'gunzip -c {gz_fpath} > {fpath}'.format(**locals())
        res = run_simple(cmdline)
        if not silent: info()
        if not res:
            return None
    return fpath

