#!/usr/bin/env python3

""" Convert Homer format to sparse upper triangular matrix """

import sys
import argparse
from collections import defaultdict
from utilities import setDefaults, createMainParent, homer2Numpy, numpy2sutm

__version__ = '1.0.0'


def homer2sutm(matrix: str, start: int, binSize: int):
    matrix = homer2Numpy(matrix)
    numpy2sutm(matrix, start, binSize)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'matrix', help='Matrix in homer format.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--start', type=int, required=True, help='Start position of matrix.')
    requiredNamed.add_argument(
        '--binSize', type=int, required=True, help='Bin size of matrix.')
    parser.set_defaults(function=homer2sutm)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
