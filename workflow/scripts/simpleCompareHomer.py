#!/usr/bin/env python3

""" Simple subtraction of 2 H5 format matrices """


import sys
import argparse
import numpy as np
from typing import List
from hicmatrix import HiCMatrix as hm
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def simpleSubtract(matrices: List, outFileName: str):

    hic1 = hm.hiCMatrix(matrices[0])
    hic2 = hm.hiCMatrix(matrices[1])

    if hic1.matrix.shape != hic2.matrix.shape:
        sys.exit("The two matrices have different size. Use matrices having "
                 "the same resolution and created using the same parameters. "
                 "Check the matrix values using the tool `hicInfo`.")

    nan_bins = set(hic1.nan_bins)
    nan_bins = nan_bins.union(hic2.nan_bins)
    new_matrix = hic1.matrix - hic2.matrix
    hic1.setMatrixValues(new_matrix)
    hic1.maskBins(sorted(nan_bins))
    hic1.save(outFileName)



def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=simpleSubtract)
    parser.add_argument('matrices', nargs=2, help='HiC matrix in homer format.')
    parser.add_argument('--outFileName', help='HiC matrix in homer format.')


    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
