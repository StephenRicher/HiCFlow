#!/usr/bin/env python3

""" Simple subtraction of 2 H5 format matrices """


import sys
import argparse
import numpy as np
from typing import List
from scipy.sparse import csr_matrix
from hicmatrix import HiCMatrix as hm
from scipy.ndimage import median_filter
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def simpleSubtract(matrices: List, outMatrix: str, outMatrixFilter: str):

    hic1 = hm.hiCMatrix(matrices[0])
    hic2 = hm.hiCMatrix(matrices[1])

    if hic1.matrix.shape != hic2.matrix.shape:
        sys.exit("The two matrices have different size. Use matrices having "
                 "the same resolution and created using the same parameters. "
                 "Check the matrix values using the tool `hicInfo`.")

    hic1.convert_to_obs_exp_matrix()
    hic2.convert_to_obs_exp_matrix()

    nan_bins = set(hic1.nan_bins)
    nan_bins = nan_bins.union(hic2.nan_bins)

    newMatrix = hic2.matrix - hic1.matrix

    for i, out in enumerate([outMatrix, outMatrixFilter]):
        if i == 1:
            filtered = median_filter(newMatrix.todense(), size=3)
            hic1.setMatrixValues(csr_matrix(filtered))
        else:
            hic1.setMatrixValues(newMatrix)
        hic1.maskBins(sorted(nan_bins))
        hic1.save(out)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=simpleSubtract)
    parser.add_argument('matrices', nargs=2, help='HiC matrix in homer format.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--outMatrix', required=True,
        help='HiC matrix in h5 format.')
    parser.add_argument(
        '--outMatrixFilter',required=True,
        help='Median filtered HiC matrix in h5 format.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
