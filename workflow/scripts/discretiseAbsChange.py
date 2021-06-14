#!/usr/bin/env python3

""" Compute absolute logFC and discretise """

import sys
import argparse
import numpy as np
import pandas as pd
from typing import List
from utilities import setDefaults, createMainParent, readHomer


__version__ = '1.0.0'


def shadowCompareHiC(matrices: List, threshold: float, maxDistance: int):

    allRegions = []
    for matrix in matrices:
        mat = readHomer(matrix, diagonal=True, sparse=False)
        chrom = mat.attrs['chrom']
        binSize = mat.attrs['binSize']

        # Remove interactions above a certain distance
        mat['seperation'] = (mat['start'] - mat['start2']).abs()
        mat = mat.loc[abs(mat['seperation']) < maxDistance]

        # Find and drop bins with a high proportion of missing values
        binZeros = mat.groupby('start')['score'].apply(
            lambda x: sum(x==0) / len(x))
        dropBins = binZeros[binZeros > threshold].index
        dropRows = mat['start'].isin(dropBins) | mat['start2'].isin(dropBins)
        mat = mat.loc[~dropRows]

        mat['abs(score)'] = abs(mat['score'])
        summed = mat.groupby('start')[['abs(score)']].sum().reset_index()

        summed['chrom'] = chrom
        summed['end'] = summed['start'] + binSize
        allRegions.append(summed)

    allRegions = pd.concat(allRegions, axis=1)
    allRegions['score'] = pd.qcut(
        allRegions['abs(score)'], 20, labels=np.linspace(0, 1, 20))

    columns = ['chrom', 'start', 'end', 'score']
    allRegions[columns].to_csv(sys.stdout, index=False, header=False, sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=shadowCompareHiC)
    parser.add_argument(
        'matrices', nargs='*', help='HiC matrix in homer format.')
    parser.add_argument(
        '--threshold', type=float, default=0.8,
        help='Mask bins with with a high proportion of zeros '
             '(default: %(default)s)')
    parser.add_argument(
        '--maxDistance', type=int, default=1000000,
        help='Remove interactions greater than this distance '
             '(default: %(default)s)')


    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
