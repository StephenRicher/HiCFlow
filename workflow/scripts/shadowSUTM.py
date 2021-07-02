#!/usr/bin/env python3

""" Randomly split and relabel HiC interaction frequencies """

import os
import sys
import argparse
import numpy as np
import pandas as pd
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def shadowSUTM(adjM: str, out: str, nShadow: int, seed: int):

    np.random.seed(seed)
    adjM = pd.read_pickle(adjM)
    allShadow = []
    for i in range(1, nShadow + 1):
        shadow = adjM.copy()
        shadow.loc[:] = shuffleAlongAxis(shadow.values, 1)
        shadow[i] = np.log2(shadow['adjIF_x'] / shadow['adjIF_y'])
        shadow = shadow.groupby('start1')[i].sum().to_frame()
        allShadow.append(shadow)
    (pd.concat(allShadow, axis=1)
        .melt(ignore_index=False, var_name='shadow', value_name='logFC')
        .sort_values('start1').to_pickle(out))


def shuffleAlongAxis(a, axis):
    """ Shuffle independenly along an axis """
    idx = np.random.rand(*a.shape).argsort(axis=axis)
    return np.take_along_axis(a, idx, axis=axis)


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
        '--seed', default=None, type=int,
        help='Seed for random swap (default: %(default)s)')
    parser.set_defaults(function=shadowSUTM)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
