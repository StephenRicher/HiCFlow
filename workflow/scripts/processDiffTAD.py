#!/usr/bin/env python3

""" Process output of differential TAD to UCSC format """

import os
import sys
import argparse
import numpy as np
import pandas as pd
from typing import List
from hicmatrix import HiCMatrix as hm
from sklearn.preprocessing import KBinsDiscretizer
from statsmodels.stats.multitest import fdrcorrection
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def processDiffTAD(
        diffTAD: List, matrix: str, outDiff: str, outPickle: str,
        adjustP: bool, alpha: float, name: str):
    allTADs = pd.concat([readDiffTAD(file) for file in diffTAD])

    pCols = ['p-value left-inter-TAD', 'p-value right-inter-TAD', 'p-value intra-TAD']
    allTADs['p'] = allTADs[pCols].fillna(1).min(axis=1)
    allTADs['p(adj)'] = fdrcorrection(allTADs['p'])[1]
    refP = 'p(adj)' if adjustP else 'p'
    allTADs['diffTAD'] = allTADs[refP] < alpha

    # Write diffTADs
    if outDiff is not None:
        cols = ['chrom', 'start', 'end']
        (allTADs.loc[allTADs['diffTAD']].sort_values(['chrom', 'start'])
            .to_csv(outDiff, sep='\t', index=False, header=False, columns=cols))


    hic = hm.hiCMatrix(matrix)
    values = np.abs(np.triu(hic.matrix.toarray()))
    allTADsummary, allSizes = processTADscores(allTADs, hic, values)
    densityStats = getDensity(allSizes, values, hic)

    mergeBy =  ['chrom', 'size']
    allTADsummary = pd.merge(densityStats, allTADsummary, left_on=mergeBy, right_on=mergeBy)
    allTADsummary['z'] = (allTADsummary['density'] - allTADsummary['mean']) / allTADsummary['std']

    mergeBy =  ['chrom', 'start', 'end']
    allTADs = pd.merge(allTADs, allTADsummary.drop('size', axis=1), left_on=mergeBy, right_on=mergeBy)

    if outPickle is not None:
        allTADs.to_pickle(outPickle)

    # Remove any domains that were too small to test
    allTADs = allTADs.loc[allTADs[['mean', 'std', 'z']].notna().all(axis=1)]
    # Bin score from 0 to 1000 for UCSC browser
    est = KBinsDiscretizer(n_bins=1000, encode='ordinal', strategy='uniform')
    allTADs['score'] = est.fit_transform(allTADs['z'].values.reshape(-1, 1))
    allTADs['strand'] = '.'
    allTADs['colour'] = 0
    allTADs['name'] = allTADs['diffTAD'].apply(lambda x: name if x else f'non-{name}')
    cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'start', 'end', 'colour']

    print(f'visibility=4 useScore="On"')
    allTADs[cols].to_csv(sys.stdout, index=False, header=False, sep='\t')


def getDensity(allSizes, values, hic):
    """ Scan all domain sizes to retrieve
        backgound mean / std of absolute change """
    chrom = hic.getChrNames()[0]
    chromSize = hic.get_chromosome_sizes()[chrom]
    binSize = hic.getBinSize()
    densityStats = {}
    for size in allSizes:
        idx1, idx2 = 0, size
        scores = []
        while True:
            scores.append(values[idx1:idx2, idx1:idx2].sum())
            idx1, idx2 = idx1 + 1, idx2 + 1
            if idx2 * binSize > chromSize:
                break
        scores = np.array(scores)
        densityStats[(chrom, size)] = (scores.mean(), scores.std())
    densityStats = pd.DataFrame(densityStats).T.reset_index()
    densityStats.columns = ['chrom', 'size', 'mean', 'std']

    return densityStats


def processTADscores(allTADs, hic, values):
    """ Loop through each TAD domains - score each domain by absolute change """
    allTADsummary = {}
    allSizes = set() # Store tad sizes
    for row in allTADs.itertuples(index=False):
        idx1, idx2 = hic.getRegionBinRange(f'{row.chrom}', row.start, row.end)
        domain = values[idx1:idx2, idx1:idx2]
        allSizes.add(len(domain)) # store matrix length
        allTADsummary[(row.chrom, row.start, row.end)] = (len(domain), domain.sum())
    allTADsummary = pd.DataFrame(allTADsummary).T.reset_index()
    allTADsummary.columns = ['chrom', 'start', 'end', 'size', 'density']

    return allTADsummary, allSizes


def readDiffTAD(file):
    dtypes = ({
        'chrom': str, 'start': int, 'end': int, 'p-value left-inter-TAD': float,
        'p-value right-inter-TAD': float, 'p-value intra-TAD': float
    })
    df = pd.read_csv(file, comment='#', usecols=[0,1,2,6,7,8],
                     names=dtypes.keys(), dtype=dtypes, sep='\t')
    return df


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=processDiffTAD)
    parser.add_argument(
        'matrix',
        help='Subtraction matrix, in .h5 format of compared matrices.')
    parser.add_argument(
        'diffTAD', nargs='+',
        help='Pair of output files from hicDifferentialTAD.')
    parser.add_argument(
        '--alpha', type=float, default=0.05,
        help='Threshold for calling diffTADs (default: %(default)s)')
    parser.add_argument(
        '--adjustP', action='store_true',
        help='If set alpha will be compared against '
             'the FDR ajusted p-value (default: %(default)s)')
    parser.add_argument(
        '--name', default='diffTAD',
        help='Name for differential TAD domains (default: %(default)s)')
    parser.add_argument(
        '--outDiff', help='Path to write BED file of diffTAD only.')
    parser.add_argument(
        '--outPickle',
        help='Path to store pickled dataframe for post processing.')


    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
