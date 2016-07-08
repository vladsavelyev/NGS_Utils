from os.path import dirname, join, abspath, splitext

from Utils.file_utils import verify_file, adjust_path
from Utils.logger import critical, debug

SUPPORTED_GENOMES = ['hg19', 'hg19-noalt', 'hg38', 'hg38-noalt', 'hg19-chr21', 'GRCh37']


def check_genome(genome):
    if genome not in SUPPORTED_GENOMES:
        critical('Genome ' + genome + ' is not supported. Supported genomes: ' + ', '.join(SUPPORTED_GENOMES))


def _get(relative_path, genome):
    # chrom = None
    # if '-chr' in genome:
    #     genome, chrom = genome.split('-')
    check_genome(genome)
    relative_path = relative_path.format(genome=genome)

    path = abspath(join(dirname(__file__), relative_path))
    # if chrom:
    #     return '<(grep {chrom} {path})'
    # else:
    return path


def get_fai(genome):
    return _get(join('fai', '{genome}.fa.fai'), genome)

def get_chrom_lengths(genome=None, fai_fpath=None):
    assert genome or fai_fpath

    if not fai_fpath:
        # if '-chr' in genome:
        #     genome, chrom = genome.split('-')
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
                    chrom, length = line.split()[0], line.split()[1]
                    chr_lengths.append((chrom, length))

    return chr_lengths

def get_chrom_order(genome=None, fai_fpath=None):
    chr_lengths = get_chrom_lengths(genome, fai_fpath)
    chr_order = {c: i for i, (c, l) in enumerate(chr_lengths)}
    return chr_order


def get_canonical_transcripts(genome):
    genome = genome.split('-')[0]
    if genome == 'GRCh37': genome = 'hg19'
    check_genome(genome)
    return _get('canonical_transcripts_{genome}.txt', genome)

def get_canonical_transcripts_ids(genome):
    with open(verify_file(get_canonical_transcripts(genome))) as f:
        return set(l.strip().split('.')[0] for l in f)
