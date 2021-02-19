#!/usr/bin/env python3

""" Convert ASHIC n*n formart to sparse upper triangular matrix """

import sys
import argparse
import itertools
import numpy as np
import pandas as pd
from collections import defaultdict
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def ashic2sutm(matrix: str, chrom: str, start: int, binSize: int,):

    # Reader homer, remove header and first 2 text rows and convert to numpy
    matrix = pd.read_csv(matrix, sep='\t', skiprows=1, header=None)
    # Must convert to integer as HiCcompare requires
    matrix = matrix.drop([0, 1], axis=1).to_numpy().astype(int)
    # Extract upper triangle only with diagonal
    matrix = np.triu(matrix, k=0)
    # Define start positions
    coords = list(range(start, (start + len(matrix) * binSize), binSize))
    # Get all binPairs
    binPairs = np.array(list(itertools.product(coords, repeat=2)))
    # Flatten to 1D list and retrieve all nonZero indexes
    matrix = matrix.flatten()
    nonZero = (matrix != 0)
    # Get bin pairs of nonZero values
    binPairs = binPairs[nonZero]
    # Remove zero values from matrix
    matrix = matrix[nonZero]
    # Write sutm to stdout
    for (bin1, bin2), value in zip(binPairs, matrix):
        print(bin1, bin2, value)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'matrix', help='Matrix in homer format.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--chrom', type=str, required=True, help='Chromosome of matrix.')
    requiredNamed.add_argument(
        '--start', type=int, required=True, help='Start position of matrix.')
    requiredNamed.add_argument(
        '--binSize', type=int, required=True, help='Bin size of matrix.')
    parser.set_defaults(function=ashic2sutm)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
