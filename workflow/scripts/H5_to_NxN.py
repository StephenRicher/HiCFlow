#!/usr/bin/env python3

""" Convert H5 matrix to NxN TSV file."""

import sys
import argparse
import numpy as np
from hicmatrix import HiCMatrix as hm
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def H5_to_NxN(matrix: str):
    matrix = hm.hiCMatrix(matrix).matrix
    np.savetxt(sys.stdout, matrix.toarray(), delimiter='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=H5_to_NxN)
    parser.add_argument('matrix', help='HiC matrix in H5 format.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
