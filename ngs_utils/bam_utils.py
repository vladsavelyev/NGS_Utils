from os.path import basename, splitext
import pysam
from ngs_utils.file_utils import verify_file, adjust_path
from ngs_utils.logger import critical, err


def verify_bam(fpath, description='', is_critical=False, silent=False):
    if not verify_file(fpath, description, is_critical=is_critical, silent=silent):
        return None

    fpath = adjust_path(fpath)

    logfn = critical if is_critical else err
    if not fpath.endswith('.bam'):
        logfn('The file ' + fpath + ' is supposed to be BAM but does not have the .bam '
            'extension. Please, make sure you pass proper file.')
        return None

    # TODO: check if binary

    return fpath


def sample_name_from_bam(in_bam):
    with pysam.Samfile(in_bam, 'rb') as bamfile:
        rg = bamfile.header.get('RG')
        if not rg or len(rg) == 0 or not rg[0].get('SM') or not rg[0].get('ID'):
            return splitext(basename(in_bam)[0])[0]
        else:
            return rg[0].get('SM', rg[0].get('ID'))
