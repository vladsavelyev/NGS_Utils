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

def get_key_genes():
    return _get('key_genes/umccr_cancer_genes.latest.tsv')

def get_key_genes_txt():
    return _get('key_genes/umccr_cancer_genes.latest.genes')

def get_key_genes_set(fpath=get_key_genes()):
    genes = set()
    with open(fpath) as f:
        for i, l in enumerate(f):
            if i != 0:
                genes.add(l.strip().split('\t')[0])
    return genes

def get_key_genes_bed(genome, is_critical=False, coding_only=False):
    return _get(f'key_genes/umccr_cancer_genes.{genome}.{"transcript" if not coding_only else "coding"}.bed',
                is_critical=is_critical)

#############################
###### Known fusions ########

"""
Bcbio uses simple-sv-annotation to prioritize SVs based on [this list of known fusions](https://github.com/AstraZeneca-NGS/simple_sv_annotation/blob/master/fusion_pairs.txt)

We extend the list with known fusions from [Hartwig's resources](https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd?path=%2FHMF-Pipeline-Resources) and feed the result into a rerun of `simple-sv-annotation`.

To make the list, first download Hartwig's fusions:

```
wget https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMF-Pipeline-Resources&files=KnownFusions.zip -O hmf.zip
unzip hmf.zip
```

Results are:

```
knownFusionPairs.csv
knownPromiscuousFive.csv
knownPromiscuousThree.csv
```

simple_sv_annotation list is coming from FusionCatcher, so we get it from there to make sure it's most recent (e.g. 11 Feb, 2019 fusioncather has 7866 pairs versus 6527 in simple_sv_annotation):

```
wget https://raw.githubusercontent.com/ndaniel/fusioncatcher/master/bin/generate_known.py
grep "        \['" generate_known.py | sed "s#        \['##" | sed "s#','#,#" | sed "s#'\],##" | sed "s#'\]##" > fusioncatcher_pairs.txt
```

Also we bring cancer genes as we have known fusion genes there:

```
cp ../key_genes/umccr_cancer_genes.latest.tsv .
```

We compare lists in `compare.R` and deside that fusionscatcher list is too big and we'll exclude it altogether.

We also compare with NGC or COSMIC Cencus fusion genes in our cancer list, they mostly overlap, so we added remaining 200 into the cancer list.
"""

def get_known_fusion_pairs():
    return _get('fusions/knownFusionPairs.csv')

def get_known_fusion_heads():
    return _get('fusions/knownPromiscuousFive.csv')

def get_known_fusion_tails():
    return _get('fusions/knownPromiscuousThree.csv')




