#!/usr/bin/env python3

""" Randomly split and relabel HiC interaction frequencies """

import os
import sys
import argparse
import numpy as np
import pandas as pd
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def shadowSUTM(mergedSUTM: str, out: str, nShadow: int, seed: int):

    np.random.seed(seed)

    mergedSUTM = pd.read_pickle(mergedSUTM)
    IF1ratio = mergedSUTM.attrs['IF1ratio']
    allShadow = []
    for i in range(nShadow):
        shadowDF = mergedSUTM.copy()
        shadowDF['adjIF_x'] = shadowDF['adjIF'].apply(
            lambda x: np.random.binomial(x, IF1ratio))
        shadowDF['adjIF_y'] = shadowDF['adjIF'] - shadowDF['adjIF_x']
        shadowDF['abs(logFC)'] = np.log2(shadowDF['adjIF_x'] / shadowDF['adjIF_y'])#.abs()
        invalid = shadowDF['abs(logFC)'].isin([np.nan, np.inf, -np.inf])
        shadowDF = shadowDF.loc[~invalid, 'abs(logFC)'].reset_index()
        shadowDF = shadowDF.groupby('start1')['abs(logFC)'].sum().reset_index()
        shadowDF['shadow'] = i + 1
        allShadow.append(shadowDF)
    pd.concat(allShadow).to_pickle(out)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'mergedSUTM', help='Pickled absolute sum logFC of merged true data.')
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
