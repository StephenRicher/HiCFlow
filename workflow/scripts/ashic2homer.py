#!/usr/bin/env python3

""" Convert ASHIC n*n formart to HOMER """

import sys
import argparse
import pandas as pd
from collections import defaultdict
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def ashic2homer(matrix: str, chrom: str, start: int, binSize: int,):

    matrix = pd.read_csv(matrix, sep='\s+', header=None).fillna(0).astype(int)
    # Define start positions
    binCoords = list(range(start, (start + len(matrix) * binSize), binSize))
    rowNames = [f'{chrom}-{pos}' for pos in binCoords]
    header = ['HiCMatrix (directory=.)', 'region'] + rowNames
    matrix.insert(0, 'rowNames1', rowNames)
    matrix.insert(0, 'rowNames2', rowNames)
    # Write homer format to stdout
    print('\t'.join(header))
    matrix.to_csv(sys.stdout, header=False, index=False, sep='\t')



def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'matrix', help='ASHIC matrix in n*n format.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--chrom', type=str, required=True, help='Chromosome of matrix.')
    requiredNamed.add_argument(
        '--start', type=int, required=True, help='Start position of matrix.')
    requiredNamed.add_argument(
        '--binSize', type=int, required=True, help='Bin size of matrix.')
    parser.set_defaults(function=ashic2homer)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
