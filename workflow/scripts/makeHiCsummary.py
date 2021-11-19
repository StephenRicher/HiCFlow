#!/usr/bin/env python3

""" Create a pseudo HiC summary file compatible with Cscore tool input """

import sys
import random
import argparse
from hicmatrix import HiCMatrix as hm
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def makeHiCsummary(matrices):
    for matrix in matrices:
        # Use path to set seed for random number generator
        random.seed(matrix)
        matrix = hm.hiCMatrix(matrix)
        binSize = matrix.getBinSize()
        start = matrix.getBinPos(0)[1] # 1-based start position
        # Ensure chrom name is prefixed with chr
        chrom = str(matrix.getChrNames()[0])
        chrom = f'chr{chrom}' if not chrom.startswith('chr') else chrom
        # Get indices of non-zero elements
        rows, columns = matrix.matrix.nonzero()
        for r, c in zip(rows, columns):
            # Don't count lower half of matrix
            if r > c:
                continue
            # Write entries according to number of contacts at that position
            for val in range(int(matrix.matrix[r, c])):
                startBase1 = start + (r * binSize)
                baseRange1 = range(startBase1, startBase1 + binSize)

                startBase2 = start + (c * binSize)
                baseRange2 = range(startBase2, startBase2 + binSize)

                # Generate pseudo-random position within bin range
                pos1 = random.choice(baseRange1)
                pos2 = random.choice(baseRange2)

                print(',', chrom, pos1, '+', chrom, pos2, '+', sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument('matrices', nargs='+', help='Matrices in h5 format.')
    parser.set_defaults(function=makeHiCsummary)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
