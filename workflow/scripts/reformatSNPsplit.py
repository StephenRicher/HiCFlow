#!/usr/bin/env python3

""" Reformat a phased VCF file to an input appropriate for SNPsplit. """

import os
import sys
import logging
import argparse
import fileinput


__version__ = '1.0.0'


def main(file, **kwargs):

    with fileinput.input(file) as fh:
        for line in fh:
            line = line.strip()
            # Skip VCF header LINES
            if line.startswith('##') or line.startswith('#CHROM'): continue
            # Extract fixed fields
            chrom, pos, id, ref, alt, qual, filter, info = line.split()[:8]
            # Extract genotype fields assuming 1 sample
            format, sample = line.split()[8:10]
            # Extract GT data
            GTfield = format.split(':').index('GT')
            GT = sample.split(':')[GTfield]
            if GT == '0|1':
                variants = f'{ref}/{alt}'
            elif GT == '1|0':
                variants = f'{alt}/{ref}'
            else:
                continue
            print(id, chrom, pos, 1, variants, sep='\t')


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument('file', help='Phased VCF file (default: stdin))')
    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'

    base = argparse.ArgumentParser(add_help=False)
    base.add_argument(
        '--version', action='version', version=f'%(prog)s {__version__}')
    base.add_argument(
        '--verbose', action='store_const', const=logging.DEBUG,
        default=logging.INFO, help='verbose logging for debugging')

    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[base, custom])
    args = parser.parse_args()

    log_format='%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
    logging.basicConfig(level=args.verbose, format=log_format)

    return args


if __name__ == '__main__':
    args = parse_arguments()
    return_code = args.function(**vars(args))
    logging.shutdown()
    sys.exit(return_code)
