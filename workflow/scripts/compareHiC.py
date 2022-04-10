#!/usr/bin/env python3

""" Simple subtraction of 2 H5 format matrices """


import sys
import argparse
import numpy as np
import pandas as pd
from typing import List
from scipy.stats import zscore
from scipy.sparse import csr_matrix
from hicmatrix import HiCMatrix as hm
from scipy.ndimage import median_filter
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def simpleSubtract(
    matrices: List, outMatrix: str, outMatrixFilter: str,
    minSum: int, raw: List):


    hic1 = hm.hiCMatrix(matrices[0])
    hic2 = hm.hiCMatrix(matrices[1])
    mask = getMask(raw, hic1, minSum)
    if (hic1.matrix.shape != hic2.matrix.shape):
        sys.exit("The two matrices have different size. Use matrices having "
                 "the same resolution and created using the same parameters. "
                 "Check the matrix values using the tool `hicInfo`.")

    nan_bins = set(hic1.nan_bins)
    nan_bins = nan_bins.union(hic2.nan_bins)
    newMatrix = (hic2.matrix - hic1.matrix).todense()

    # On occasion the raw may have excess rows - remove these
    mask = mask[np.indices(newMatrix.shape, sparse=True)]

    for i, out in enumerate([outMatrixFilter, outMatrix]):
        if i == 0:
            filtered = median_filter(newMatrix, size=3)
            filtered[mask] = 0
            hic1.setMatrixValues(csr_matrix(filtered))
        else:
            newMatrix[mask] = 0
            hic1.setMatrixValues(newMatrix)
        hic1.maskBins(sorted(nan_bins))
        hic1.save(out)


def getMask(raw, hic, minSum=0):
    if raw:
        raw1 = hm.hiCMatrix(raw[0])
        raw2 = hm.hiCMatrix(raw[1])
        mask = (raw1.matrix + raw2.matrix).todense() < minSum
    else:
        mask = hic.matrix.todense() < 0
    return mask


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=simpleSubtract)
    parser.add_argument('matrices', nargs=2, help='HiC matrix in homer format.')
    parser.add_argument(
        '--minSum', type=int, default=0,
        help='Total per cell of raw matrices must be atleast '
             'this value (default: %(default)s).')
    parser.add_argument(
        '--raw', nargs=2, help='Raw matrices to filter by minimum bin count.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--outMatrix', required=True, help='HiC matrix in h5 format.')
    requiredNamed.add_argument(
        '--outMatrixFilter',required=True,
        help='Median filtered HiC matrix in h5 format.')


    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
