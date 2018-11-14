from os.path import dirname, join, abspath, splitext, isfile

from ngs_utils.file_utils import verify_file, adjust_path, verify_dir
from ngs_utils.key_genes_utils import get_genes_from_file
from ngs_utils.logger import critical, debug


SUPPORTED_GENOMES = ['hg19', 'hg19-noalt', 'hg38', 'hg38-noalt', 'hg19-chr21', 'GRCh37', 'mm10']


def check_genome(genome):
    if genome not in SUPPORTED_GENOMES:
        critical('Genome ' + str(genome) + ' is not supported. Supported genomes: ' + ', '.join(SUPPORTED_GENOMES))

def _get(relative_path, genome=None, is_critical=False):
    if genome:
        check_genome(genome)
    else:
        genome = ''
    relative_path = relative_path.format(genome=genome)

    path = abspath(join(dirname(__file__), relative_path))
    if is_critical:
        return verify_file(path, is_critical=True)
    return path


######################
######## FAI #########
def get_fai(genome, is_critical=False):
    return _get(join('fai', '{genome}.fa.fai'), genome, is_critical=is_critical)

def get_chrom_lengths(genome=None, fai_fpath=None):
    assert genome or fai_fpath, f'One of genome or fai_fpath should be not None: genome={genome}, fai_fpath={fai_fpath}'

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

def ucsc_to_ensembl(genome, is_critical=False):
    """ mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" hg19 > hg19.ucscToEnsembl.tsv
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" hg38 > hg38.ucscToEnsembl.tsv
    """
    return _get(join('fai', '{genome}.ucscToEnsembl.tsv'), genome, is_critical=is_critical)


#############################
######## Gene lists #########

def get_signatures_probabilities(is_critical=False):
    return _get('signatures_probabilities.txt', is_critical=is_critical)

def get_suppressors(is_critical=False):
    return _get('suppressors.txt', is_critical=is_critical)

def get_key_genes_bed(genome, is_critical=False, coding_only=False):
    return _get(f'key_genes/key_genes.{genome}.{"transcript" if not coding_only else "coding"}.bed',
                is_critical=is_critical)

def get_key_genes_set(is_critical=False):
    return get_genes_from_file(_get('key_genes/key_genes.txt', is_critical=is_critical))
