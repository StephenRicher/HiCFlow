#!/usr/bin/env python3

""" Calculate adjusted absolute logFC from adjIF1 and adjIF2 SUTM matrices
    computed by HiCcompare. Also calculate per-bin merged counts for procesing
    by shadowSUTM.py. """

import os
import sys
import pickle
import argparse
import numpy as np
import pandas as pd
from typing import List
from matplotlib import cm
from scipy.stats import wilcoxon
from matplotlib.colors import to_hex
from statsmodels.stats.multitest import fdrcorrection
from utilities import setDefaults, createMainParent
try:
    from pandarallel import pandarallel
except ModuleNotFoundError:
    pass


__version__ = '1.0.0'


def computeAdjM(
        adjIF1: str, adjIF2: str, nShadow: int, fdr: float, cmap: str,
        rawOut: str, binSize: int, chrom: str, seed: int, threads: int):

    np.random.seed(seed)

    adjIF1 = readSUTM(adjIF1, lower=True)
    adjIF2 = readSUTM(adjIF2, lower=True)

    adjM = pd.merge(adjIF1, adjIF2, how='inner', left_index=True, right_index=True)
    # Removal uninformative bins where IF is equal
    adjM = adjM.loc[adjM['adjIF_x'] != adjM['adjIF_y']]

    # Perform Wilcoxon between adjIF
    direction = adjM.groupby('start1').apply(directionPreference).rename(
        {0: 'p_Wilcoxon', 1: 'direction'}, axis=1)

    # Get per-bin direction and set to -1 or 1
    adjM['diff'] = (adjM['adjIF_x'] - adjM['adjIF_y']).apply(getDirection)

    if (threads > 1) and ('pandarallel' in sys.modules):
        pandarallel.initialize(nb_workers=5, verbose=0)
        permuted = adjM.groupby('start1').parallel_apply(permuteTest, nShadow)
    else:
        permuted = adjM.groupby('start1').apply(permuteTest, nShadow)
    permuted = permuted.rename(
        {0: 'p_Permute', 1: 'score', 2: 'length'}, axis=1)

    # Merge permutation stats with Wilcoxon
    permuted = pd.merge(permuted, direction, left_index=True, right_index=True).reset_index()

    permuted['p(adj)_Wilcoxon'] = fdrcorrection(permuted['p_Wilcoxon'])[1]
    permuted['p(adj)_Permute'] = fdrcorrection(permuted['p_Permute'])[1]
    permuted['pScore'] = -np.log10(permuted['p_Permute'])

    # Compute minimum p-value of permutation test
    maxScore = -np.log10(1 / nShadow)
    permuted['colour'] = permuted.apply(
        getColour, args=(maxScore, fdr, cmap), axis=1)
    permuted['chrom'] = chrom
    permuted['end'] = permuted['start1'] + binSize
    if rawOut:
        permuted.to_csv(rawOut, index=False, header=False, sep='\t')

    # Write direction score as a BED file
    permuted['name'] = '.'
    permuted['strand'] = '.'
    permuted['thickStart'] = permuted['start1']
    permuted['thickEnd'] = permuted['end']
    columns = ([
        'chrom', 'start1', 'end', 'name', 'pScore', 'strand',
        'thickStart', 'thickEnd', 'colour'])
    # Remove sequences of length 1 (score = na)
    permuted = permuted.loc[permuted['length'] > 1]
    permuted[columns].to_csv(sys.stdout, index=False, header=False, sep='\t')


def getColour(x, maxScore, fdr, colourmap):
    """ Define RGB colours based on significance of Wilcoxon result.
        Non-significant bins are coloured on a binary scale according
        to their pScore. Bins with significance directional preference
        are also coloured according to direction. """
    # Scale pScores between 0 and 1 for indexing cmap
    pScore_scale = scale01(x['pScore'], 0, maxScore)
    if x['p(adj)_Wilcoxon'] <= fdr:
        if x['direction'] == 1:
            i = (pScore_scale * 0.5)  + 0.5
        else:
            i = (1 - pScore_scale) * 0.5
        colour = to_hex(cm.get_cmap(colourmap, 100)(i))[1:]
    else:
        i = pScore_scale
        colour = to_hex(cm.get_cmap('binary', 100)(i))[1:]
    colour = f'{int(colour[:2], 16)},{int(colour[2:4], 16)},{int(colour[4:], 16)}'
    return colour


def permuteTest(x, nShadow):
    values = x['diff'].values
    score = float(patternSum(np.array([values])))
    boolRand = (np.random.random((nShadow, len(x))) > 0.5).astype(int)
    randRLE = patternSum(boolRand)
    totalAbove = (randRLE >= score).sum()
    p = (totalAbove / nShadow)
    if p == 0:
        p += (1 / nShadow)
    return pd.Series([p, score, len(values)])


def directionPreference(x):
    _, p = wilcoxon(x['adjIF_x'], x['adjIF_y'], alternative='two-sided')
    # Find which median is greater to determine direction
    direction = x['adjIF_x'].median() < x['adjIF_y'].median()
    direction = 1 if direction else -1
    return pd.Series([p, direction])


def patternSum(x, axis=1):
    """ Given a sequence of 1s and 0s compute the cumulative
        sum of its differences. Same as previous scores +1
        different scores -1 """
    # Find positions of change (e.g. +1 or -1)
    diff = np.diff(x, axis)
    # Score no change as 1, change as -1
    scores = np.where(diff != 0, -1, 1)
    sumScores = scores.sum(axis)
    # Rescale between 0 and 1
    size = np.size(x, axis)
    maxScore = (size - 1) # e.g. all same sequence
    minScore = -maxScore # e.g. alternating sequence
    return scale01(sumScores, minScore, maxScore)


def scale01(x, minX, maxX):
    """ Scale x between 0 and 1 """
    return (x - minX) / (maxX - minX)


def readSUTM(sutm, lower=False):
    sutm = pd.read_csv(sutm, names=['start1', 'start2', 'adjIF'], sep=' ')
    if lower:
        sltm = sutm.loc[sutm['start1'] != sutm['start2']].rename(
            {'start1': 'start2', 'start2': 'start1'}, axis=1)
        sutm  = pd.concat([sutm, sltm])
    return sutm.sort_values(['start1', 'start2']).set_index(['start1', 'start2'])


def getDirection(x):
    """ Return sign of number, 0s are not expected """
    return -1 if x < 0 else 1


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'adjIF1', help='Sparse upper triangular matrix file, group 1.')
    parser.add_argument(
        'adjIF2', help='Sparse upper triangular matrix file, group 2.')
    parser.add_argument(
        '--rawOut', help='Output file for unfiltered results')
    parser.add_argument(
        '--nShadow', default=1000, type=int,
        help='Number of random permutations (default: %(default)s)')
    parser.add_argument(
        '--fdr', default=0.05, type=float,
        help='FDR threshold for determing Wilcoxon '
             'significance (default: %(default)s)')
    parser.add_argument(
        '--cmap', default='bwr',
        help='Matplotlib divering colourmap for colouring '
             'Wilcoxon bed track (default: %(default)s)')
    parser.add_argument(
        '--chrom', help='Reference of correponding matrices.')
    parser.add_argument(
        '--binSize', type=int, help='Bin size of corresponding matrices')
    parser.add_argument(
        '--threads', default=1, type=int,
        help='Threads for parallel processing (default: %(default)s)')
    parser.add_argument(
        '--seed', type=int,
        help='Seed for random number generator (default: %(default)s)')
    parser.set_defaults(function=computeAdjM)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
