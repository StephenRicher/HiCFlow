#!/usr/bin/env python3

""" Generate equal length BED file of windowed region
    size for input to CscoreTool. """


import os
import sys
import argparse
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def generateEqualLengthBed(regions: str, window: int):
    with open(regions) as fh:
        for line in fh:
            chrom, start, end = line.strip().split()[:3]
            start = int(start)
            end = int(end)
            size = end - start
            nWindows = int((size + window - 1) / window)
            for i in range(nWindows):
                lower = i * window + 1
                upper = (i + 1) * window
                if upper > size:
                    upper = size
                print(chrom, lower, upper, sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument('regions', metavar='REGIONs',
        help='Regions to split in BED format.')
    parser.add_argument('window', metavar='WINDOW-SIZE', type=int,
        help='Length of split BED intervals.')
    parser.set_defaults(function=generateEqualLengthBed)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
