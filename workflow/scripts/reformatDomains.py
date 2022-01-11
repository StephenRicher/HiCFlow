#!/usr/bin/env python3

""" Reformat OnTAD output BED file to links scaled to correct coordinates."""

import sys
import logging
import argparse
import fileinput
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def reformatDomains(bed: str):

    with fileinput.input(bed) as fh:
        for line in fh:
            # Skip header line
            if line.startswith(('#', 'track')):
                continue
            chrom, start, end = line.split()[:3]
            print(chrom, start, end, '.', '0', '.', start, end, '.', sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument(
        'bed', nargs='?', default=[],
        help='BED format TAD files for conversion (default: stdin)')
    parser.set_defaults(function=reformatDomains)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
