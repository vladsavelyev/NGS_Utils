#!/usr/bin/env python
import sys
from optparse import OptionParser, SUPPRESS_HELP
from os.path import dirname, join

import GeneAnnotation
from Utils.annotate_bed import annotate
from Utils.bam_bed_utils import verify_bed
from Utils.file_utils import verify_file
from Utils.logger import critical, info


def main():
    if len(sys.argv[1]) < 0:
        critical('Usage: ' + __file__ + ' Input_BED_file -g hg19 -o Annotated_BED_file')
    input_bed_fpath = verify_bed(sys.argv[1], is_critical=True, description='Input BED file for ' + __file__)

    options = [
        (['-o', '--output-file'], dict(
            dest='output_file',
            metavar='FILE',
            help='Output file',
        )),
        (['--reference', '--features'], dict(
            dest='features',
        )),
        (['--reuse'], dict(
            dest='reuse_intermediate',
            help='reuse intermediate non-empty files in the work dir from previous run',
            action='store_true',
        )),
        (['--work-dir'], dict(dest='work_dir', metavar='DIR', help=SUPPRESS_HELP)),
        (['--log-dir'], dict(dest='log_dir', metavar='DIR', help=SUPPRESS_HELP)),
        (['-g', '--genome'], dict(dest='genome')),
    ]

    parser = OptionParser(description='Annotating BED file based on reference features annotations.')
    for args, kwargs in options:
        parser.add_option(*args, **kwargs)
    opts, args = parser.parse_args()

    features_fpath = None
    if opts.features:
        features_fpath = verify_bed(opts.features, is_critical=True)
    elif opts.genome:
        features_fpath = verify_file(join(dirname(GeneAnnotation.__file__), GeneAnnotation.ANNOTATION.format(genome=opts.genome)))
        if not features_fpath:
            critical('Genome ' + opts.genome + ' is not supported. Supported: ' + ', '.join(GeneAnnotation.SUPPORTED_GENOMES))
    else:
        critical('Error: neither --features nor --genome is specified. Features are required for annotation.')

    annotate(input_bed_fpath, features_fpath, opts.output_file)

    info('Done.')


if __name__ == '__main__':
    main()
