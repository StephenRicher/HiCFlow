#!/usr/bin/env python3

""" Compute Z score of logFC for each HiC interval by generating
   randomised matrices """

import sys
import random
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
from utilities import setDefaults, createMainParent, readHomer


__version__ = '1.0.0'


def shadowCompareHiC(file: str, nShadow: int, seed: int):

    random.seed(seed)

    mat = readHomer(file, diagonal=True, sparse=False)
    chrom = mat.attrs['chrom']
    binSize = mat.attrs['binSize']

    mat['seperation'] = (mat['start'] - mat['start2']).abs()
    mat['abs(score)'] = abs(mat['score'])
    summed = mat.groupby('start')[['abs(score)']].sum()

    # Compute mean and standard deviation per randomised bin in loop
    summed['shadowSum'] = 0
    summed['shadowSumX2'] = 0
    for _ in range(nShadow):
        mat['randomScore'] = (
            mat.groupby('seperation')['abs(score)']
               .transform(np.random.permutation))
        shadowMat = mat.groupby('start')['randomScore'].sum()
        summed['shadowSum'] += shadowMat
        summed['shadowSumX2'] += shadowMat * shadowMat
    summed = summed.reset_index()

    summed['shadowMean'] = summed['shadowSum'] / nShadow
    summed['shadowStd'] = (np.sqrt((summed['shadowSumX2'] / nShadow)
                           - (summed['shadowMean'] * summed['shadowMean'])))
    summed['Z'] = (summed['abs(score)'] - summed['shadowMean']) / summed['shadowStd']
    summed['p'] = stats.norm.sf(summed['Z'])
    summed['p(adj)'] = fdrcorrection(summed['p'])[1]
    summed['pScore'] = -10 * np.log10(summed['p(adj)'])
    summed['chrom'] = chrom
    summed['end'] = summed['start'] + binSize

    columns = ['chrom', 'start', 'end', 'pScore']
    summed[columns].to_csv(sys.stdout, index=False, header=False, sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=shadowCompareHiC)
    parser.add_argument('file', help='HiC matrix in homer format.')
    parser.add_argument(
        '--nShadow', type=int, default=100,
        help='Number of shuffled matrices to generate '
             'Z distribution (default: %(default)s)')
    parser.add_argument(
        '--seed', type=int, default=None,
        help='Seed for random number generator (default: %(default)s)')


    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
