#!/usr/bin/env python3

""" Randomly split and relabel HiC interaction frequencies """

import os
import sys
import pickle
import argparse
import numpy as np
import pandas as pd
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def shadowSUTM(adjM: str, out: str, nShadow: int, nSplits: int, seed: int):

    np.random.seed(seed)
    adjM = pd.read_pickle(adjM)
    allShadow = []
    for i in range(1, nShadow + 1):
        shadow = adjM.copy()
        shadow.loc[:] = shuffleAlongAxis(shadow.values, 1)
        shadow[i] =  shadow['adjIF_x'] - shadow['adjIF_y']
        # np.log2(shadow['adjIF_x'] / shadow['adjIF_y'])
        shadow = shadow.groupby('start1')[i].sum().abs().to_frame()
        allShadow.append(shadow)
    allShadow = (pd.concat(allShadow, axis=1)
        .melt(ignore_index=False, var_name='shadow', value_name='logFC')
        .sort_values('start1'))

    # Write data in groups sorted by start positions. This allows chunksize
    # processing of data during permution test to save memory.
    allShadow['split'] = splitData(len(allShadow), nShadow, nSplits)
    with open(out,'wb') as fh:
        for split, splitGroup in allShadow.groupby('split'):
            pickle.dump(splitGroup.drop('split', axis=1), fh)


def shuffleAlongAxis(a, axis):
    """ Shuffle independenly along an axis """
    idx = np.random.rand(*a.shape).argsort(axis=axis)
    return np.take_along_axis(a, idx, axis=axis)


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


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'adjM', help='Pickled merged IF of HiC groups.')
    parser.add_argument(
        '--out', help='Pickled output of shadow logFC.')
    parser.add_argument(
        '--nShadow', default=1, type=int,
        help='Number of random permutations (default: %(default)s)')
    parser.add_argument(
        '--nSplits', default=1, type=int,
        help='Number of splits to save pickled data (default: %(default)s)')
    parser.add_argument(
        '--seed', default=None, type=int,
        help='Seed for random swap (default: %(default)s)')
    parser.set_defaults(function=shadowSUTM)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
