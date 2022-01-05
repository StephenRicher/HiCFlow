#!/usr/bin/env python3

""" Convert H5 matrix to NxN3p TSV file."""

import sys
import argparse
import numpy as np
from hicmatrix import HiCMatrix as hm
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def H5_to_NxN3p(matrix):
    hic = hm.hiCMatrix(matrix)
    chrom = hic.getChrNames()[0]
    binSize = hic.getBinSize()
    matrix = hic.matrix.toarray()
    # Looping over numpy is slow but more memory efficient than
    # converting to pandas to append columns
    for i, row in enumerate(matrix):
        start = i * binSize
        end = start + binSize
        row = '\t'.join(row.astype(str))
        print(f'{chrom}\t{start}\t{end}\t{row}')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=H5_to_NxN3p)
    parser.add_argument('matrix', help='HiC matrix in H5 format.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
