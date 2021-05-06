#!/usr/bin/env python3

""" Pad an NxN matrix to encompass all bins in a reference. Used
    by HiCFlow to reformat matrices for matrix2Butlr script. """

import sys
import math
import logging
import argparse
import numpy as np
import pandas as pd
from utilities import setDefaults, createMainParent, readChromSizes

__version__ = '1.0.0'


def padNxN(matrix: str, ref: str, start: int, binSize: int, refSizes: str):

    # Read matrix 1 row at a time
    matrix = pd.read_csv(matrix, header=None, sep='\t', dtype=str, chunksize=1)
    # Retrieve first row of NxN matrix
    row = matrix.get_chunk()
    # Number of bins in unpadded matrix
    nValidBins = row.shape[1]
    # Get range of defined matrix bins
    startBin = start // binSize
    endBin = startBin + nValidBins
    # Get reference size in base pairs for matrix
    refSize = readChromSizes(refSizes)[ref]
    # List of all bins for padded matrix
    allBins = range(0, math.ceil(refSize / binSize))
    emptyRow = '\t'.join(['0'] * len(allBins))
    # Create empty bins to pad either side of region
    nPadLeft = startBin
    leftPad = '\t'.join(['0'] * nPadLeft) + '\t' if nPadLeft > 0 else ''
    nPadRight = allBins[-1] - endBin + 1
    rightPad = '\t'.join(['0'] * nPadRight)
    for bin1 in allBins:
        # Whole row empty
        if (bin1 < startBin) or (bin1 >= endBin):
            print(emptyRow)
        else:
            countValues = '\t'.join(row.values.tolist()[0])
            # Add a tab if there is right padding
            countValues = countValues + '\t' if nPadRight > 0 else countValues
            print(leftPad, countValues, rightPad, sep='')
            try:
                row = matrix.get_chunk()
            except StopIteration:
                pass


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument('matrix', help='HiC matrix in NxN format.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--ref', type=str, required=True,
        help='Reference / chromsome of matrix')
    requiredNamed.add_argument(
        '--start', type=int, required=True, help='Start position of matrix.')
    requiredNamed.add_argument(
        '--binSize', type=int, required=True,
        help='Bin size of matrix in base pairs.')
    requiredNamed.add_argument(
        '--refSizes', type=str, required=True,
        help='Reference / chromsome sizes file.')
    parser.set_defaults(function=padNxN)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
