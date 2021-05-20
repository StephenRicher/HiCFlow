#!/usr/bin/env python3

""" Convert TAD domain BED file to TAD boundary BED  """

import os
import sys
import math
import logging
import argparse
import fileinput
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def TADdomains2boundaries(binSize: int, bed: str):
    halfBin = int(binSize / 2)
    with fileinput.input(bed) as fh:
        for line in fh:
            chrom, start, end = line.split()[:3]
            start = int(start)
            end = int(end)
            print(chrom, start - halfBin, start + halfBin, sep='\t')
            print(chrom, end - halfBin, end + halfBin, sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument(
        'binSize', type=int,
        help='Boundary size, centered at the boundary of a TAD.')
    parser.add_argument('bed', nargs='?', default=[], help='Input BED file.')
    parser.set_defaults(function=TADdomains2boundaries)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
