#!/usr/bin/env python3

""" Process output of differential TAD to UCSC format """

import os
import sys
import logging
import argparse
import numpy as np
import pandas as pd
from typing import List
from hicmatrix import HiCMatrix as hm
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def processDiffTAD(
        diffTAD: List, subtractionMat: str, countMats: List, outDiff: str,
        outPickle: str, threshold: float, name: str):

    allTADs = pd.concat([readBED(file) for file in diffTAD])

    hic = hm.hiCMatrix(subtractionMat)
    chrom = hic.getChrNames()[0]

    # Only keep domains matching the provided matrix
    allTADs = allTADs.loc[allTADs['chrom'] == chrom]

    # Convert genomic intervals to list of tuples
    allTADs = list(allTADs.to_records(index=False))

    values = np.abs(np.triu(hic.matrix.toarray()))
    allTADsummary, allSizes = processTADscores(allTADs, hic, values)
    backgroundStats = getDomainBackground(allSizes, values, hic, countMats)

    mergeBy =  ['chrom', 'size']
    allTADs = pd.merge(
        backgroundStats, allTADsummary, left_on=mergeBy, right_on=mergeBy)
    allTADs['z'] = (allTADs['sumDiff'] - allTADs['mean']) / allTADs['std']

    # Remove any domains that were too small to test
    allTADs = allTADs.loc[allTADs[['mean', 'std', 'z']].notna().all(axis=1)]

    allTADs['rank'] = allTADs['z'].rank(pct=True)

    # Define diffTAD as top 'threshold' % of domains
    allTADs['diffTAD'] = (allTADs['rank'] * 100) >= 100 - threshold

    # Write diffTADs
    if outDiff is not None:
        cols = ['chrom', 'start', 'end']
        (allTADs.loc[allTADs['diffTAD']].sort_values(['chrom', 'start'])
            .to_csv(outDiff, sep='\t', index=False, header=False, columns=cols))

    if outPickle is not None:
        allTADs.to_pickle(outPickle)

    # Set score between 0 - 1000 according to rank
    allTADs['colour'] = 0
    allTADs['strand'] = '.'
    allTADs['score'] = allTADs['rank'] * 1000
    allTADs['name'] = allTADs['diffTAD'].apply(
        lambda x: name if x else f'non-{name}')
    cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'start', 'end', 'colour']

    print(f'visibility=4 useScore="On"')
    allTADs[cols].to_csv(sys.stdout, index=False, header=False, sep='\t')


def getDomainBackground(allSizes, values, hic, countMats):
    """ Scan all domain sizes to retrieve
        backgound mean / std of absolute change """
    chrom = hic.getChrNames()[0]
    chromSize = hic.get_chromosome_sizes()[chrom]
    binSize = hic.getBinSize()
    densityStats = {}

    IF1 = hm.hiCMatrix(countMats[0]).matrix.toarray()
    IF2 = hm.hiCMatrix(countMats[1]).matrix.toarray()
    sparsity = (IF1 + IF2) > 0 # False at empty

    for size in allSizes:
        idx1, idx2 = 0, size
        scores = []
        while True:
            score = values[idx1:idx2, idx1:idx2].sum()
            isZero = sparsity[idx1:idx2, idx1:idx2].sum()
            idx1, idx2 = idx1 + 1, idx2 + 1
            if True: #isZero > 0:
                scores.append(score)
            if idx2 * binSize > chromSize:
                break
        scores = np.array(scores)
        densityStats[(f'{chrom}', size)] = (scores.mean(), scores.std())
    densityStats = pd.DataFrame(densityStats).T.reset_index()
    densityStats.columns = ['chrom', 'size', 'mean', 'std']

    return densityStats


def processTADscores(allTADs, hic, values):
    """ Loop through each TAD domains - score each domain by absolute change """
    allTADsummary = {}
    allSizes = set() # Store tad sizes
    for chrom, start, end in allTADs:
        try:
            idx1, idx2 = hic.getRegionBinRange(f'{chrom}', start, end)
        except TypeError:
            logging.error(f'Skipping {chrom}:{start}-{end} - out of range.')
            continue
        domain = values[idx1:idx2, idx1:idx2]
        allSizes.add(len(domain)) # store matrix length
        allTADsummary[(chrom, start, end)] = (len(domain), domain.sum())
    allTADsummary = pd.DataFrame(allTADsummary).T.reset_index()

    allTADsummary.columns = ['chrom', 'start', 'end', 'size', 'sumDiff']

    return allTADsummary, allSizes


def readBED(file):
    dtypes = {'chrom': str, 'start': int, 'end': int}
    df = pd.read_csv(file, comment='#', usecols=[0,1,2],
                     names=dtypes.keys(), dtype=dtypes, sep='\t')
    return df


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=processDiffTAD)
    parser.add_argument(
        'subtractionMat',
        help='Subtraction matrix, in .h5 format of compared matrices.')
    parser.add_argument(
        'countMats', nargs=2,
        help='Count matrices of compared samples, used to exclude '
             'completely empty regions')
    parser.add_argument(
        'diffTAD', nargs='+',
        help='Pair of output files from hicDifferentialTAD.')
    parser.add_argument(
        '--threshold', type=float, default=10,
        help='Define top % of domains as differential (default: %(default)s)')
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
