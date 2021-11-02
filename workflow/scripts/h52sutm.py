#!/usr/bin/env python3

""" Convert H5 format to sparse upper triangular matrix """

import sys
import argparse
import pandas as pd
from hicmatrix import HiCMatrix as hm
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def h52sutm(matrix):
    hic = hm.hiCMatrix(matrix)
    nonzeroIdx = hic.matrix.nonzero()
    nonzeroValues = hic.matrix[nonzeroIdx].tolist()[0]
    mat = pd.DataFrame({
        'start' : nonzeroIdx[0],
        'start2': nonzeroIdx[1],
        'value' : nonzeroValues
    })
    mat = mat.loc[mat['start'] >= mat['start2']]
    mat.to_csv(sys.stdout, header=False, index=False, sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument('matrix', help='Matrix in H5 format.')
    parser.set_defaults(function=h52sutm)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
