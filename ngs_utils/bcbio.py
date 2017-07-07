import re
import six
import yaml
import sys
from os import listdir
from os.path import join, abspath, pardir, splitext, basename, dirname, realpath, isdir, isfile

from ngs_utils.Sample import BaseSample
from ngs_utils.bam_utils import verify_bam
from ngs_utils.call_process import run
from ngs_utils.config import load_yaml_config
from ngs_utils.file_utils import adjust_path, verify_dir, file_exists, safe_mkdir, verify_file, add_suffix
from ngs_utils.logger import critical, debug, info, err, warn

from ngs_reporting.utils import get_target_genes, is_small_target

import variant_filtering as vf


CALLER = 'vardict'


class BcbioSample(BaseSample):
    def __init__(self):
        BaseSample.__init__(self)
        self.old_name = None
        self.raw_name = None
        self.bcbio_project = None
        self.project_tag = None
        self.genome_build = None
        self.sample_info = None
        self.fastq_files = None

        self.is_rnaseq = None
        self.is_wgs = None
        self.min_allele_fraction = None
        self.coverage_interval = None

        self.variant_regions_bed = None
        self.sv_regions_bed = None
        self.coverage_bed = None

        self.bam = None
        self.counts_file = None

        self.batch = None
        self.batch_names = []
        self.phenotype = None

    def get_name_for_files(self):  # In case if the sample if symlink from another project, and the name was changed in this one
        return self.old_name or self.name

    def load_from_sample_info(self, sample_info, bcbio_project, need_coverage_interval=True, need_vardict=True):
        self.sample_info = sample_info
        self.bcbio_project = bcbio_project
        self.raw_name = str(sample_info['description'])
        self.name = self.raw_name.replace('.', '_')
        if 'description_original' in sample_info:
            self.old_name = str(sample_info['description_original']).replace('.', '_')
        self.fastq_files = sample_info['files']
        self.dirpath = verify_dir(join(self.bcbio_project.final_dir, self.name))
        if not verify_dir(self.dirpath):
            critical('Directory for sample ' + self.name + ' is not found in ' + self.bcbio_project.final_dir)
        self.var_dirpath = join(self.dirpath, BcbioProject.var_dir)

        self.genome_build = sample_info['genome_build']
        self.variant_regions_bed = self.bcbio_project.config_path(val=sample_info['algorithm'].get('variant_regions'))
        self.sv_regions_bed = self.bcbio_project.config_path(val=sample_info['algorithm'].get('sv_regions')) or self.variant_regions_bed
        self.coverage_bed = self.bcbio_project.config_path(val=sample_info['algorithm'].get('coverage')) or self.sv_regions_bed
        self.is_rnaseq = 'rna' in sample_info['analysis'].lower()
        self.min_allele_fraction = (1.0/100) * float(sample_info['algorithm'].get('min_allele_fraction', 1.0))
        # if not self.is_rnaseq and need_coverage_interval:
            # self.coverage_interval = _parse_coveragre_interval(self.bcbio_project, self.name)
            # if self.coverage_interval is None:
            #     warn('Not found covearge interval for sample ' + self.name)
            # if self.coverage_interval not in ['genome', 'regional', 'amplicon']:
            #     err('Error: coverage_interval is ' + str(self.coverage_interval))
        if self.variant_regions_bed is None:
            self.coverage_interval = 'genome'
        else:
            self.coverage_interval = 'regional'
        self.is_wgs = self.coverage_interval == 'genome'

        bam = adjust_path(join(self.dirpath, self.get_name_for_files() + '-ready.bam'))
        if isfile(bam) and verify_bam(bam):
            self.bam = bam
        else:
            warn('No BAM file found for ' + self.name)

        if self.is_rnaseq:
            gene_counts = adjust_path(join(self.dirpath, self.get_name_for_files() + '-ready.counts'))
            if isfile(gene_counts) and verify_file(gene_counts):
                self.counts_file = gene_counts
            else:
                warn('Counts for ' + self.name + ' not found')

        else:
            self.phenotype = sample_info.get('metadata', dict()).get('phenotype') or 'tumor'
            batch_info = sample_info.get('metadata', dict()).get('batch') or self.get_name_for_files() + '-batch'
            self.batch_names = batch_info.split(', ') if isinstance(batch_info, six.string_types) else batch_info
            if len(self.batch_names) > 1 and self.phenotype != 'normal':
                critical('Multiple batches for non-normal ' + self.phenotype + ' sample ' + self.name + ': ' + ', '.join(self.batch_names))

            variantcallers = sample_info['algorithm'].get('variantcaller') or []
            if isinstance(variantcallers, six.string_types):
                variantcallers = [variantcallers]
            if CALLER not in variantcallers:
                warn('Warning: "' + CALLER + '" is not in the variant callers (' + str(variantcallers) + ')')

    def find_mutation_files(self, passed=True):
        return _find_mutation_files(join(self.dirpath, BcbioProject.varfilter_dir), passed)

    def find_raw_vcf(self, silent=False):
        vcf_fpath = None
        if self.batch and self.phenotype != 'normal':
            vcf_fpath = self.bcbio_project.find_vcf_file(self.batch.name, silent=silent)
        if not vcf_fpath:  # in sample dir?
            if not silent:
                debug('-')
                debug('Not found VCF in the datestamp dir, looking at the sample-level dir')
                debug('-')
            vcf_fpath = self.bcbio_project.find_vcf_file_from_sample_dir(self, silent=silent or self.phenotype == 'normal')
        return vcf_fpath
        
    def find_annotated_vcf(self):
        return verify_file(join(self.dirpath, BcbioProject.varannotate_dir,
                                self.name + '-' + CALLER + BcbioProject.anno_vcf_ending + '.gz'), silent=True)

    def find_filt_vcf(self, passed=False):
        path = join(self.dirpath, BcbioProject.varfilter_dir, self.name + '-' + CALLER +
                    ((BcbioProject.filt_vcf_ending + '.gz') if not passed else BcbioProject.pass_filt_vcf_ending))
        return verify_file(path, silent=True)

    def find_mutation_file(self, passed=True):
        mut_fname = CALLER + '.' + vf.mut_file_ext
        mut_fpath = join(self.dirpath, BcbioProject.varfilter_dir, mut_fname)
        if passed:
            mut_fpath = add_suffix(mut_fpath, vf.mut_pass_suffix)
        return verify_file(mut_fpath, silent=True)

    def find_sv_vcf(self):
        return self.find_cnv_file(self.name + '-manta.vcf.gz') or \
               self.find_cnv_file(self.name + '-lumpy.vcf.gz')

    def find_sv_tsv(self):
        return self.find_cnv_file(self.name + '-sv-prioritize.tsv')

    def find_seq2c_calls(self):
        return self.find_cnv_file(self.name + '-seq2c.tsv') or \
              (self.find_cnv_file(self.batch.name + '-seq2c.tsv') if self.batch else None) or \
               verify_file(join(self.bcbio_project.date_dir, BcbioProject.cnv_dir, BcbioProject.seq2c_fname), silent=True) or \
               verify_file(join(self.bcbio_project.date_dir, BcbioProject.cnv_dir, 'Seq2C.tsv'), silent=True)

    def find_seq2c_coverage(self):
        return self.find_cnv_file(self.name + '-seq2c-coverage.tsv') or \
              (self.find_cnv_file(self.batch.name + '-seq2c-coverage.tsv') if self.batch else None) or \
               verify_file(join(self.bcbio_project.date_dir, BcbioProject.cnv_dir, 'seq2c-cov.tsv'), silent=True) or \
               verify_file(join(self.bcbio_project.date_dir, BcbioProject.cnv_dir, 'cov.tsv'), silent=True)

    def find_cnvkit_file(self):
        return self.find_cnv_file(self.name + '-cnvkit.cnr') or \
              (self.find_cnv_file(self.batch.name + '-cnvkit.cnr') if self.batch else None)

    def find_seq2c_file(self):
        return self.find_cnv_file(self.name + '-seq2c.tsv') or \
              (self.find_cnv_file(self.batch.name + '-seq2c.tsv') if self.batch else None)

    def find_cnv_file(self, fname):
        for fpath in [join(self.dirpath, fname),
                      join(self.dirpath, BcbioProject.cnv_dir, fname)]:
            if isfile(fpath):
                return verify_file(fpath, silent=True)

    def find_coverage_stats(self):
        return verify_file(join(self.dirpath, 'qc', 'coverage', self.name + '_coverage.bed'), silent=True)
    
    def find_ngs_report(self, silent=False):
        return \
            verify_file(join(self.bcbio_project.date_dir, BcbioProject.reports_dir,
                             self.name + '.html')) or \
            verify_file(join(self.dirpath, BcbioProject.ngs_report_name,
                             BcbioProject.ngs_report_name + '.html'), silent=silent)

    def get_metric(self, names):
        if isinstance(names, six.string_types):
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
        return self.get_metric('Total_reads')
    
    def is_dedupped(self):
        return self.sample_info.get('algorithm', {}).get('mark_duplicates', False)
    
    
class Batch:
    def __init__(self, name=None):
        self.name = name
        self.normal = None
        self.tumor = []

    def __str__(self):
        return self.name


class BcbioProject:
    varfilter_dir = 'varFilter'
    varannotate_dir = 'varAnnotate'
    cnv_dir = 'cnv'
    var_dir = 'var'
    ngs_report_name = 'ngs_report'
    reports_dir = 'reports'
    anno_vcf_ending = '.anno.vcf'
    filt_vcf_ending = '.anno.filt.vcf'
    pass_filt_vcf_ending = '.anno.filt.' + vf.mut_pass_suffix + '.vcf'
    seq2c_fname = 'seq2c.tsv'
    cnvkit_fname = 'cnvkit.tsv'
    evaluate_panel_dir = 'eval_panel'
    oncoprints_dir = 'oncoprints'

    ## RNAseq
    counts_names = ['counts.tsv', 'dexseq.tsv', 'gene.sf.tpm.tsv', 'isoform.sf.tpm.tsv']
    expression_dir = 'expression'

    multiqc_report_name = 'report.html'
    call_vis_name = 'call_vis.html'

    def __init__(self):
        self.dir = None
        self.config_dir = None
        self.final_dir = None
        self.date_dir = None
        self.log_dir = None
        self.postproc_log_dir = None
        self.work_dir = None

        self.var_dir = None
        self.raw_var_dir = None
        self.expression_dir = None
        self.raw_expression_dir = None

        self.samples = []
        self.batch_by_name = dict()

        self.variant_regions_bed = None
        self.sv_regions_bed = None          # "sv_regions" or "variant_regions"
        self.original_coverage_bed = None
        self.coverage_bed = None            # "coverage" or "sv_regions" or "variant_regions"

        self.project_name = None
        self.genome_build = None
        self.coverage_interval = None  # amplicon, regional, genome
        self.min_allele_fraction = None
        self.is_wgs = None
        self.is_rnaseq = None

    def set_project_level_dirs(self, bcbio_cnf, project_name=None, final_dir=None, create_dirs=False,
                               proc_name='postproc'):
        assert self.dir

        self.final_dir = self.set_final_dir(bcbio_cnf, self.dir, final_dir)
        if create_dirs: safe_mkdir(self.final_dir)

        self.project_name = self._set_project_name(self.dir, project_name)

        self.work_dir = abspath(join(self.dir, 'work', proc_name))
        if create_dirs: safe_mkdir(self.work_dir)

        self.date_dir = self._set_date_dir(bcbio_cnf, self.final_dir, create_dir=create_dirs)
        self.log_dir = join(self.date_dir, 'log')
        self.postproc_log_dir = join(self.log_dir, proc_name)
        if create_dirs: safe_mkdir(self.postproc_log_dir)

        self.var_dir = join(self.date_dir, BcbioProject.var_dir)
        self.raw_var_dir = join(self.var_dir, 'raw')
        self.expression_dir = join(self.date_dir, BcbioProject.expression_dir)
        self.raw_expression_dir = join(self.expression_dir, 'raw')

    def load_from_bcbio_dir(self, input_dir, project_name=None, proc_name='postproc',
                            need_coverage_interval=True, need_vardict=True):
        """
        Analyses existing bcbio folder.
        :param input_dir: root bcbio folder, or any other directory inside it
        """
        self.dir, detected_final_dir = detect_bcbio_dir(input_dir)
        self.config_dir = abspath(join(self.dir, 'config'))
        bcbio_cnf, _ = load_bcbio_cnf(self.config_dir)
        self.set_project_level_dirs(bcbio_cnf, project_name=project_name, final_dir=detected_final_dir,
                                    proc_name=proc_name)
        self.set_samples(bcbio_cnf, need_coverage_interval=need_coverage_interval, need_vardict=need_vardict)
        self._load_bcbio_summary()
        # self._load_target_info()

    def set_samples(self, bcbio_cnf, need_coverage_interval=True, need_vardict=True):
        debug('Reading sample details...')
        for sample_info in bcbio_cnf['details']:
            s = BcbioSample()
            s.load_from_sample_info(sample_info, bcbio_project=self,
                                    need_coverage_interval=need_coverage_interval, need_vardict=need_vardict)
            self.samples.append(s)
        if any(not s.bam for s in self.samples):
            warn('ERROR: for some samples, BAM files not found in the final dir')

        self.samples.sort(key=lambda _s: _s.key_to_sort())
        self.batch_by_name = self.update_batches(self.samples)

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
            info('RNAseq')
        elif self.coverage_interval:
            info('Coverage interval: ' + str(self.coverage_interval))

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
                s.sample_info['metrics'] = metrics_by_sample[s.name]
            
    # def _load_target_info(self):
    #     interval = None
    #     target_info_file = join(self.date_dir, 'multiqc', 'report', 'metrics', 'target_info.yaml')
    #     if not isfile(target_info_file):
    #         target_info_file = join(self.date_dir, 'log', 'multiqc_bcbio', 'report', 'metrics', 'target_info.yaml')
    #     if isfile(target_info_file):
    #         debug('Parsing ' + target_info_file + ' for coverage_interval')
    #         with open(target_info_file) as f:
    #             d = yaml.load(f)
    #             if 'coverage_interval' in d:
    #                 interval = d['coverage_interval']
    #             else:
    #                 err('Not found coverage_interval in ' + target_info_file)
    #     return interval
    
    # def init_new(self, proj_dir, bcbio_cnf, project_name=None):
    #     """ Creates folders
    #     """
    #     self.dir = proj_dir
    #     self.config_dir = safe_mkdir(abspath(join(self.dir, 'config')))
    #     with open(join(self.config_dir, 'bcbio.yaml'), 'w') as yaml_file:
    #         yaml_file.write(save_yaml(bcbio_cnf))
    #
    #     self._set_props(bcbio_cnf, project_name=project_name, create_dirs=True)
    #     self._set_samples(bcbio_cnf)

    def config_path(self, val):
        if not val:
            return val
        return adjust_path(join(self.config_dir, val))

    @staticmethod
    def _set_date_dir(bcbio_cnf, final_dirpath, create_dir=False):
        if 'fc_date' not in bcbio_cnf:
            critical('Error: fc_date not in bcbio config!')
        if 'fc_name' not in bcbio_cnf:
            critical('Error: fc_name not in bcbio config!')
        # Date dirpath is from bcbio and named after fc_name, not our own project name
        date_dirpath = join(final_dirpath, bcbio_cnf['fc_date'] + '_' + bcbio_cnf['fc_name'])
        if create_dir: safe_mkdir(date_dirpath)
        elif not verify_dir(date_dirpath):
            critical('Error: no project directory of format {fc_date}_{fc_name}')
        return date_dirpath

    @staticmethod
    def _set_project_name(dirpath, project_name=None):
        small_project_path = None
        if '/ngs/oncology/analysis/' in realpath(dirpath):
            short_path = realpath(dirpath).split('/ngs/oncology/analysis/')[1]  # bioscience/Bio_0031_Heme_MRL_DLBCL_IRAK4/bcbio_Dev_0079
            small_project_path = '/'.join(short_path.split('/')[1:])
        if not project_name:
            # path is like /ngs/oncology/analysis/bioscience/Bio_0031_Heme_MRL_DLBCL_IRAK4/bcbio_Dev_0079
            if small_project_path:
                project_name = '_'.join(small_project_path.split('/'))  # Bio_0031_Heme_MRL_DLBCL_IRAK4_bcbio_Dev_0079
        second_part = basename(dirpath)  # bcbio_Dev_0079
        bcbio_project_parent_dirname = basename(dirname(dirpath))  # Bio_0031_Heme_MRL_DLBCL_IRAK4
        if not project_name:
            project_name = bcbio_project_parent_dirname + '_' + second_part
        return project_name

    @staticmethod
    def set_final_dir(bcbio_cnf, proj_dir, final_dir=None, create_dir=False):
        if final_dir:
            return final_dir
        elif 'upload' in bcbio_cnf and 'dir' in bcbio_cnf['upload']:
            final_dirname = bcbio_cnf['upload']['dir']
            final_dir = adjust_path(join(proj_dir, 'config', final_dirname))
            if create_dir: safe_mkdir(final_dir)
            verify_dir(final_dir, 'upload directory specified in the bcbio config', is_critical=True)
        else:
            final_dir = join(proj_dir, 'final')
            if create_dir: safe_mkdir(final_dir)
            if not verify_dir(final_dir):
                critical('If final directory it is not named "final", please, specify it in the bcbio config.')
        return final_dir
    
    @staticmethod
    def update_batches(samples):
        batch_by_name = {bn: Batch(bn) for bn in list(set([b for s in samples for b in s.batch_names]))}
        for sample in samples:
            for bn in sample.batch_names:
                batch_by_name[bn].name = bn
                sample.batch = batch_by_name[bn]
                if sample.phenotype == 'normal':
                    if batch_by_name[bn].normal:
                        critical('Multiple normal samples for batch ' + bn)
                    batch_by_name[bn].normal = sample
                else:
                    batch_by_name[bn].tumor.append(sample)

        for batch in batch_by_name.values():
            if batch.normal and not batch.tumor:
                info('Batch ' + batch.name + ' contains only normal, treating sample ' + batch.normal.name + ' as tumor')
                batch.normal.phenotype = 'tumor'
                batch.normal.batch = batch
                batch.tumor = [batch.normal]
                batch.normal = None

        # setting up batch properties
        for b in batch_by_name.values():
            if b.normal and b.tumor:
                b.paired = True
            else:
                b.paired = False
        for b in batch_by_name.values():
            for t_sample in b.tumor:
                t_sample.normal_match = b.normal

        return batch_by_name

    def find_vcf_file(self, batch_name, silent=False):
        vcf_fname = batch_name + '-' + CALLER + '.vcf'
        annot_vcf_fname = batch_name + '-' + CALLER + '-annotated.vcf'

        vcf_annot_fpath_gz = adjust_path(join(self.date_dir, annot_vcf_fname + '.gz'))  # in datestamp
        var_raw_vcf_annot_fpath_gz = adjust_path(join(self.raw_var_dir, annot_vcf_fname + '.gz'))  # in datestamp/var/raw

        vcf_fpath_gz = adjust_path(join(self.date_dir, vcf_fname + '.gz'))  # in datestamp
        var_vcf_fpath_gz = adjust_path(join(self.var_dir, vcf_fname + '.gz'))  # in datestamp/var
        var_raw_vcf_fpath_gz = adjust_path(join(self.raw_var_dir, vcf_fname + '.gz'))  # in datestamp/var/raw

        vcf_fpath = adjust_path(join(self.date_dir, vcf_fname))  # in datestamp
        var_vcf_fpath = adjust_path(join(self.var_dir, vcf_fname))  # in datestamp/var
        var_raw_vcf_fpath = adjust_path(join(self.raw_var_dir, vcf_fname))  # in datestamp/var/raw

        if isfile(vcf_annot_fpath_gz):
            verify_file(vcf_annot_fpath_gz, is_critical=True)
            info('Found annotated VCF in the datestamp dir ' + vcf_annot_fpath_gz)
            return vcf_annot_fpath_gz
        else:
            debug('Not found annotated VCF in the datestamp dir ' + vcf_annot_fpath_gz)

        if isfile(var_raw_vcf_annot_fpath_gz):
            verify_file(var_raw_vcf_annot_fpath_gz, is_critical=True)
            info('Found annotated VCF in the datestamp/var/raw dir ' + var_raw_vcf_annot_fpath_gz)
            return var_raw_vcf_annot_fpath_gz
        else:
            debug('Not found annotated VCF in the datestamp/var/raw dir ' + var_raw_vcf_annot_fpath_gz)

        if isfile(vcf_fpath_gz):
            verify_file(vcf_fpath_gz, is_critical=True)
            info('Found VCF in the datestamp dir ' + vcf_fpath_gz)
            return vcf_fpath_gz
        else:
            debug('Not found VCF in the datestamp dir ' + vcf_fpath_gz)

        if isfile(var_raw_vcf_fpath_gz):
            verify_file(var_raw_vcf_fpath_gz, is_critical=True)
            info('Found VCF in the datestamp/var/raw dir ' + var_raw_vcf_fpath_gz)
            return var_raw_vcf_fpath_gz
        else:
            debug('Not found VCF in the datestamp/var/raw dir ' + var_raw_vcf_fpath_gz)

        if isfile(vcf_fpath):
            verify_file(vcf_fpath, is_critical=True)
            info('Found uncompressed VCF in the datestamp dir ' + vcf_fpath)
            return vcf_fpath
        else:
            debug('Not found uncompressed VCF in the datestamp dir ' + vcf_fpath)

        if isfile(var_raw_vcf_fpath):
            verify_file(var_raw_vcf_fpath, is_critical=True)
            info('Found uncompressed VCF in the datestamp/var/raw dir ' + var_raw_vcf_fpath)
            return var_raw_vcf_fpath
        else:
            debug('Not found uncompressed VCF in the datestamp/var/raw dir ' + var_raw_vcf_fpath)

        if isfile(var_vcf_fpath_gz):
            verify_file(var_vcf_fpath_gz, is_critical=True)
            info('Found VCF in the datestamp/var dir ' + var_vcf_fpath_gz)
            return var_vcf_fpath_gz
        else:
            debug('Not found VCF in the datestamp/var dir ' + var_vcf_fpath_gz)

        if isfile(var_vcf_fpath):
            verify_file(var_vcf_fpath, is_critical=True)
            info('Found uncompressed VCF in the datestamp/var dir ' + var_vcf_fpath)
            return var_vcf_fpath
        else:
            debug('Not found uncompressed VCF in the datestamp/var dir ' + var_vcf_fpath)

        if not silent:
            warn('Warning: no VCF found for batch ' + batch_name + ', ' + CALLER + ', gzip or '
                'uncompressed version in the datestamp directory.')
        return None

    @staticmethod
    def find_vcf_file_from_sample_dir(sample, silent=False):
        vcf_fname = sample.get_name_for_files() + '-' + CALLER + '.vcf'
        
        sample_var_dirpath = join(sample.dirpath, 'var')
        vcf_fpath_gz = adjust_path(join(sample.dirpath, vcf_fname + '.gz'))  # in var
        var_vcf_fpath_gz = adjust_path(join(sample_var_dirpath, vcf_fname + '.gz'))  # in var
        var_raw_vcf_fpath_gz = adjust_path(join(sample_var_dirpath, 'raw', vcf_fname + '.gz'))  # in var
        vcf_fpath = adjust_path(join(sample.dirpath, vcf_fname))
        var_vcf_fpath = adjust_path(join(sample_var_dirpath, vcf_fname))  # in var
        var_raw_vcf_fpath = adjust_path(join(sample_var_dirpath, 'raw', vcf_fname))  # in var

        if isfile(vcf_fpath_gz):
            verify_file(vcf_fpath_gz, is_critical=True)
            info('Found VCF ' + vcf_fpath_gz)
            return vcf_fpath_gz
        else:
            debug('Not found VCF ' + vcf_fpath_gz)

        if isfile(var_vcf_fpath_gz):
            verify_file(var_vcf_fpath_gz, is_critical=True)
            info('Found VCF in the var/ dir ' + var_vcf_fpath_gz)
            return var_vcf_fpath_gz
        else:
            debug('Not found VCF in the var/ dir ' + var_vcf_fpath_gz)

        if isfile(var_raw_vcf_fpath_gz):
            verify_file(var_raw_vcf_fpath_gz, is_critical=True)
            info('Found VCF in the var/raw/ dir ' + var_raw_vcf_fpath_gz)
            return var_raw_vcf_fpath_gz
        else:
            debug('Not found VCF in the var/raw/ dir ' + var_raw_vcf_fpath_gz)

        if isfile(vcf_fpath):
            verify_file(vcf_fpath, is_critical=True)
            info('Found uncompressed VCF ' + vcf_fpath)
            return vcf_fpath
        else:
            debug('Not found uncompressed VCF ' + vcf_fpath)

        if isfile(var_vcf_fpath):
            verify_file(var_vcf_fpath, is_critical=True)
            info('Found uncompressed VCF in the var/ dir ' + var_vcf_fpath)
            return var_vcf_fpath
        else:
            debug('Not found VCF in the var/ dir ' + var_vcf_fpath)

        if isfile(var_raw_vcf_fpath):
            verify_file(var_raw_vcf_fpath, is_critical=True)
            info('Found uncompressed VCF in the var/raw/ dir ' + var_raw_vcf_fpath)
            return var_raw_vcf_fpath
        else:
            debug('Not found VCF in the var/raw/ dir ' + var_raw_vcf_fpath)

        if not silent:
            warn('Warning: no VCF found for ' + sample.name + ', ' + CALLER + ', gzip or uncompressed version in and outside '
                'the var directory. Phenotype is ' + str(sample.phenotype))
        return None

    def find_seq2c_file(self):
        return verify_file(join(self.date_dir, BcbioProject.cnv_dir, BcbioProject.seq2c_fname), silent=True) or \
               verify_file(join(self.date_dir, BcbioProject.cnv_dir, 'Seq2C.tsv'), silent=True)

    def find_seq2c_filt_file(self):
        return verify_file(join(self.date_dir, BcbioProject.cnv_dir,
                                add_suffix(BcbioProject.seq2c_fname, 'filt')), silent=True)

    def find_seq2c_coverage(self):
        return verify_file(join(self.date_dir, BcbioProject.cnv_dir, 'seq2c-cov.tsv'), silent=True)

    def find_cnvkit_file(self):
        return verify_file(join(self.date_dir, BcbioProject.cnv_dir, BcbioProject.cnvkit_fname), silent=True)

    def find_cnvkit_filt_file(self):
        return verify_file(join(self.date_dir, BcbioProject.cnv_dir,
                                add_suffix(BcbioProject.cnvkit_fname, 'filt')), silent=True)

    def get_cnv_caller(self):
        seq2c_file = self.find_seq2c_filt_file()
        cnvkit_file = self.find_cnvkit_filt_file()
        if not seq2c_file or (cnvkit_file and self.is_wgs):
            return 'CNVkit'
        else:
            return 'Seq2C'

    def find_cnv_filt_file(self):
        caller = self.get_cnv_caller()
        if caller == 'Seq2C':
            return self.find_seq2c_filt_file()
        else:
            return self.find_cnvkit_filt_file()

    def find_multiqc_report(self):
        for fpath in [
            join(self.date_dir, BcbioProject.multiqc_report_name),
            join(self.date_dir, 'multiqc_postproc', 'multiqc_report.html'),
        ]:
            if verify_file(fpath, silent=True):
                return fpath

    def find_mutation_files(self, passed=True):
        return _find_mutation_files(self.var_dir, passed=passed)
    
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

    def get_target_genes(self):
        return get_target_genes(self.genome_build, self.coverage_bed)
        
    def is_small_target(self):
        return is_small_target(self.coverage_bed)


def _find_mutation_files(base_dir, passed=True):
    mut_fname = CALLER + '.' + vf.mut_file_ext
    mut_fpath = join(base_dir, mut_fname)
    single_mut_fpath = add_suffix(mut_fpath, vf.mut_single_suffix)
    paired_mut_fpath = add_suffix(mut_fpath, vf.mut_paired_suffix)
    fpaths = [mut_fpath, single_mut_fpath, paired_mut_fpath]
    if passed:
        fpaths = [add_suffix(p, vf.mut_pass_suffix) for p in fpaths]
    return [p for p in fpaths if verify_file(p, silent=True)]


def detect_bcbio_dir(input_dir, silent=False):
    """
    :param input_dir: root dir, or any other directory inside it
    :return: (root_dir, final_dir=None)
    """
    root_dir, final_dir = None, None
    input_dir = abspath(input_dir)

    if isdir(join(input_dir, 'config')):
        root_dir = input_dir
        final_dir = None

    elif isdir(abspath(join(input_dir, pardir, 'config'))):
        root_dir = abspath(join(input_dir, pardir))
        if 'final' in basename(input_dir):
            final_dir = input_dir

    elif isdir(abspath(join(input_dir, pardir, pardir, 'config'))):
        root_dir = abspath(join(input_dir, pardir, pardir))
        if 'final' in basename(abspath(join(input_dir, pardir))):
            final_dir = abspath(join(input_dir, pardir))

    elif isdir(abspath(join(input_dir, pardir, pardir, pardir, 'config'))):
        root_dir = abspath(join(input_dir, pardir, pardir, pardir))
        if 'final' in basename(abspath(join(input_dir, pardir, pardir))):
            final_dir = abspath(join(input_dir, pardir, pardir))
    else:
        if silent:
            return None, None
        else:
            critical(
                'Are you running on a bcbio directory?\n'
                'Can\'t find `config` directory at ' + join(input_dir, 'config') + ' or ' + abspath(join(input_dir, pardir, 'config')) + '. '
                'Make sure that you changed to a bcbio root or final directory, or provided it as a first argument.\n'
                'Type `' + basename(__file__) + ' --help` for more help.')
    
    if not silent:
        info('Bcbio project directory: ' + root_dir)
        if final_dir:
            info('"final" directory: ' + final_dir)
    return root_dir, final_dir


def load_bcbio_cnf(config_dir):
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
    info('Using bcbio YAML config ' + yaml_fpath)
    return load_yaml_config(yaml_fpath), yaml_fpath


def _normalize(name):
    return name.lower().replace('_', '').replace('-', '')


def ungzip_if_needed(cnf, fpath):
    if fpath.endswith('.gz'):
        fpath = fpath[:-3]
    if not file_exists(fpath) and file_exists(fpath + '.gz'):
        gz_fpath = fpath + '.gz'
        cmdline = 'gunzip -c {gz_fpath}'.format(**locals())
        res = run(cmdline, output_fpath=fpath)
        info()
        if not res:
            return None
    return fpath

 
# def _parse_coveragre_interval(bcbio_project, sample_name):
#     interval = None
#     target_info_file = join(bcbio_project.date_dir, 'multiqc', 'report', 'metrics', 'target_info.yaml')
#     if not isfile(target_info_file):
#         target_info_file = join(bcbio_project.date_dir, 'log', 'multiqc_bcbio', 'report', 'metrics', 'target_info.yaml')
#     if isfile(target_info_file):
#         debug('Parsing ' + target_info_file + ' for coverage_interval')
#         with open(target_info_file) as f:
#             d = yaml.load(f)
#             if 'coverage_interval' in d:
#                 interval = d['coverage_interval']
#             else:
#                 debug('Not found coverage_interval in ' + target_info_file)
#
#     if not interval:
#         bcbio_log = bcbio_project.find_in_log('bcbio-nextgen.log', silent=False)
#         if not bcbio_log:
#             debug('Not found bcbio-nextgen.log')
#         else:
#             debug('Parsing coverage_interval from log')
#             pattern = r'.* (?P<sample>.+): Assigned coverage as \'(?P<interval>.+)\''
#             with open(bcbio_log) as f:
#                 for l in f:
#                     m = re.match(pattern, l)
#                     if m and len(m.groups()) == 2:
#                         if m.group('sample') == sample_name:
#                             interval = m.group('interval')
#             if not interval:
#                 err('Coverage interval info is not found in ' + bcbio_log)
#     return interval
