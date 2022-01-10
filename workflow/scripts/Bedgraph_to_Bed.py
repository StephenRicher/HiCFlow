#!/usr/bin/env python3

""" Convert Bedgraph to Bed """

import sys
import argparse
import fileinput
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'

def Bedgraph_to_Bed(bedgraph):
    with fileinput.input(bedgraph) as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('track '):
                continue
            chrom, start, end, score = line.split()
            print(chrom, start, end, ".", score, sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=Bedgraph_to_Bed)
    parser.add_argument('bedgraph', nargs='?', default=[],
        help='Input Bedgraph file (default: stdin)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
