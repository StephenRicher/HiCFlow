#!/usr/bin/env python3

""" Reformat a phased VCF file to an input appropriate for SNPsplit. """

import sys
import logging
import argparse
import fileinput
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def reformatSNPsplit(VCF: str, flipSNP: bool):

    with fileinput.input(VCF) as fh:
        for line in fh:
            line = line.strip()
            # Skip VCF header
            if line.startswith('##') or line.startswith('#CHROM'):
                continue
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
            if flipSNP:
                a1, a2 = variants.split('/')
                variants = f'{a2}/{a1}'
            print(id, chrom, pos, 1, variants, sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=reformatSNPsplit)
    parser.add_argument(
        'VCF', nargs='?', help='Phased VCF file (default: stdin))')
    parser.add_argument(
        '--flipSNP', action='store_true',
        help='Flip allelic phased SNP assignment, this is '
             'often arbitrary (default: %(default)s)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
