#!/usr/bin/env python3

""" Reformat OnTAD output BED file to links scaled to correct coordinates."""

import sys
import logging
import argparse
import fileinput
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def reformatDomains(bed: str, scale: int, trimChr: bool, **kwargs):

    with fileinput.input(bed) as fh:
        for line in fh:
            # Skip header line
            if line.startswith(('#', 'track')):
                continue
            columns = line.split()
            chrom = columns[0]
            if chrom.startswith('chr') and trimChr:
                chrom = chrom[3:]
            start = str(int(columns[1]) + scale)
            end = str(int(columns[2]) + scale)

            print(chrom, start, end, '.', '0', '.', start, end, '.', sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument(
        'bed', nargs='?', default=[],
        help='BED format TAD files for conversion (default: stdin)')
    parser.add_argument(
        '--scale', type=int, default=0,
        help='0-based genomic coordinate of start position '
             'to rescale (default: %(default)s)')
    parser.add_argument(
        '--trimChr', action='store_true',
        help='Trim chr prefix from chrom name (default: %(default)s)')
    parser.set_defaults(function=reformatDomains)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
