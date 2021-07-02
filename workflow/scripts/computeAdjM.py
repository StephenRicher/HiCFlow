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
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def computeAdjM(adjIF1: str, adjIF2: str, adjM_out: str, merge_out: str, nSplits: int):

    adjIF1 = readSUTM(adjIF1, lower=True)
    adjIF2 = readSUTM(adjIF2, lower=True)

    # Compute absolute logFC
    adjM = pd.merge(
        adjIF1, adjIF2, how='inner',
        left_index=True, right_index=True)
    # Write merged data to pickle for fast access later
    adjM.to_pickle(merge_out)


    adjM[0] = (adjM['adjIF_x'] - adjM['adjIF_y']).apply(getDirection)

    # Sum values across bins and save to pickle
    adjM = (adjM.groupby('start1')[0].apply(rle).to_frame()
        .melt(ignore_index=False, var_name='shadow', value_name='logFC')
        .sort_values('start1'))

    # Write data in groups sorted by start positions. This allows chunksize
    # processing of data during permution test to save memory. Must match
    # the number of splits used in shadowSUTM.py
    adjM['split'] = splitData(len(adjM), 1, nSplits)
    with open(adjM_out,'wb') as fh:
        for split, splitGroup in adjM.groupby('split'):
            pickle.dump(splitGroup.drop('split', axis=1), fh)


def splitData(size, nShadow, nSplits):
    """ Split start positions as evenly as possible"""
    uniqueStart = int(size / nShadow)
    # Cannot split more than unique start as shadows must be grouped
    assert nSplits <= uniqueStart
    if nSplits == 1:
        return np.repeat(1, size)
    else:
        splitSize = uniqueStart // nSplits
        evenSplits = np.repeat(np.arange(1, nSplits), splitSize)
        # Create final split with remainder
        remaining = uniqueStart - (splitSize * (nSplits - 1))
        finalSplits = np.ones(remaining) * nSplits
        # Merge even with remainder
        splits = np.append(evenSplits, finalSplits)
    return np.repeat(splits, nShadow)


def readSUTM(sutm, lower=False):
    sutm = pd.read_csv(sutm, names=['start1', 'start2', 'adjIF'], sep=' ')
    if lower:
        sltm = sutm.loc[sutm['start1'] != sutm['start2']].rename(
            {'start1': 'start2', 'start2': 'start1'}, axis=1)
        sutm  = pd.concat([sutm, sltm])
    return sutm.sort_values(['start1', 'start2']).set_index(['start1', 'start2'])


def getDirection(x):
    if x < 0:
        return -1
    elif x > 0:
        return 1
    else:
        return 0


def rle(inarray):
    """ run length encoding. Partial credit to R rle function.
        Multi datatype arrays catered for including non Numpy
        returns: tuple (runlengths, startpositions, values) """
    # https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encodi
    ia = np.asarray(inarray)
    n = len(ia)
    if n == 0:
        return None
    else:
        y = ia[1:] != ia[:-1]               # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)   # must include last element posi
        z = np.diff(np.append(-1, i))       # run lengths
        return len(z)


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
        '--adjM_out',
        help='Pickled output for absolute sum logFC.')
    parser.add_argument(
        '--merge_out', help='Pickled output of merged per-bin counts.')
    parser.add_argument(
        '--nSplits', default=1, type=int,
        help='Number of splits to save pickled data (default: %(default)s)')
    parser.set_defaults(function=computeAdjM)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
