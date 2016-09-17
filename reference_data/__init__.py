from genericpath import isfile
from os.path import dirname, join, abspath, splitext

from Utils.file_utils import verify_file, adjust_path, verify_dir
from Utils.logger import critical, debug

SUPPORTED_GENOMES = ['hg19', 'hg19-noalt', 'hg38', 'hg38-noalt', 'hg19-chr21', 'GRCh37', 'mm10']

def check_genome(genome):
    if genome not in SUPPORTED_GENOMES:
        critical('Genome ' + str(genome) + ' is not supported. Supported genomes: ' + ', '.join(SUPPORTED_GENOMES))

def _get(relative_path, genome=None):
    if genome:
        check_genome(genome)
    else:
        genome = ''
    relative_path = relative_path.format(genome=genome)

    path = abspath(join(dirname(__file__), relative_path))
    return path

######################
### TRICKY_REGIONS ###
def find_tricky_regions_dir(genome):
    short_genome = genome.split('-')[0]
    dirpath = _get(join('tricky_regions', '{genome}'), short_genome)
    if not verify_dir(dirpath):
        not_found_err('directory ' + dirpath)
    return dirpath

def not_found_err(what):
    critical('Tricky regions ' + what + ' not found. To resolve this, please ' +
             'either install git-lfs (https://git-lfs.github.com) and pull the repository again, ' +
             'or provide the path to the tricky regions explicitly using the --tricky-regions option')

def find_tricky_regions(dirpath, name):
    return verify_file(join(dirpath, name + '.bed.gz'), is_critical=True)

tricky_regions_full_names = {
    'bad_promoter': 'Bad promoter',
    'gc0to15': 'GC 0-15%',
    'gc15to20': 'Low GC 15-20%',
    'gc20to25': 'Low GC 20-25%',
    'gc25to30': 'Low GC 25-30%',
    'gc65to70': 'High GC 65-70%',
    'gc70to75': 'High GC 70-75%',
    'gc75to80': 'High GC 75-80%',
    'gc80to85': 'High GC 80-85%',
    'gc85to100': 'High GC 85-100%',
    'low_complexity_lt51bp': 'Low complexity <51bp',
    'low_complexity_51to200bp': 'Low complexity 51-200bp',
    'low_complexity_gt200bp': 'Low complexity >200bp',
    'heng_universal_mask': 'Heng\'s mask',
    'repeats': 'Repeats',
    'self_chain': 'Self chain',
}

######################
######## FAI #########
def get_fai(genome):
    return _get(join('fai', '{genome}.fa.fai'), genome)

def get_chrom_lengths(genome=None, fai_fpath=None):
    assert genome or fai_fpath, 'One of genome or fai_fpath should be not None: ' \
                                'genome=' + str(genome) + ' fai_fpath=' + str(fai_fpath)

    if not fai_fpath:
        check_genome(genome)
        fai_fpath = get_fai(genome)
    else:
        fai_fpath = verify_file(fai_fpath, is_critical=True)
        if not fai_fpath.endswith('.fai') and not fai_fpath.endswith('.fa'):
            critical('Error: .fai or .fa is accepted.')

    chr_lengths = []

    if fai_fpath.endswith('.fa'):
        debug('Reading genome sequence (.fa) to get chromosome lengths')
        with open(fai_fpath, 'r') as handle:
            from Bio import SeqIO
            reference_records = SeqIO.parse(handle, 'fasta')
            for record in reference_records:
                chrom = record.id
                chr_lengths.append((chrom, len(record.seq)))

    else:
        debug('Reading genome index file (.fai) to get chromosome lengths')
        with open(fai_fpath, 'r') as handle:
            for line in handle:
                line = line.strip()
                if line:
                    chrom, length = line.split()[0], int(line.split()[1])
                    chr_lengths.append((chrom, length))

    return chr_lengths

def get_chrom_order(genome=None, fai_fpath=None):
    chr_lengths = get_chrom_lengths(genome, fai_fpath)
    chr_order = {c: i for i, (c, l) in enumerate(chr_lengths)}
    return chr_order

#########################
###### TRANSCRIPTS ######
def get_canonical_transcripts(genome):
    short_genome = genome.split('-')[0]
    if short_genome == 'GRCh37':
        short_genome = 'hg19'
    check_genome(short_genome)
    return _get(join('canonical_transcripts', 'canonical_transcripts_{genome}.txt'), short_genome)

def get_canonical_transcripts_ids(genome):
    fpath = get_canonical_transcripts(genome)
    if not isfile(fpath):
        return None
    with open(verify_file(fpath, is_critical=True)) as f:
        return set(l.strip().split('.')[0] for l in f)

###################
###### GENES ######
def key_genes_800():
    return _get('az_key_genes.800.txt')
