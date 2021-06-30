#!/usr/bin/env python3

""" Perform permutation testing by merging BAMs from both groups and randomly
    split (keeping the relative sample sizes approximately equal). Output is
    written to SUTM format and is repeated a determined number of times. """

import os
import sys
import argparse
import numpy as np
import pandas as pd
from typing import List
from collections import defaultdict
from utilities import setDefaults, createMainParent, getGroup

__version__ = '1.0.0'


def homer2sutm(SUTM1: str, SUTM2: str, out1: str, out2: str, seed: int):

    np.random.seed(seed)

    SUTM1 = pd.read_csv(SUTM1, names=['pos1', 'pos2', 'count'], sep='\s+')
    SUTM2 = pd.read_csv(SUTM2, names=['pos1', 'pos2', 'count'], sep='\s+')

    totalHiC = SUTM1['count'].sum() + SUTM2['count'].sum()
    SUTM1ratio = SUTM1['count'].sum() / totalHiC

    mergedSUTM = (
        pd.concat([SUTM1, SUTM2]).groupby(['pos1', 'pos2']).sum().reset_index())

    mergedSUTM['SUTM1'] = mergedSUTM['count'].apply(
        lambda x: np.random.binomial(x, SUTM1ratio))
    mergedSUTM['SUTM2'] = mergedSUTM['count'] - mergedSUTM['SUTM1']

    for group in ['SUTM1', 'SUTM2']:
        names = ['pos1', 'pos2', group]
        out = out1 if group == 'SUTM1' else out2
        sutm = mergedSUTM.loc[mergedSUTM[group] > 0, names].sort_values(['pos1', 'pos2'])
        sutm.to_csv(out, sep=' ', header=False, index=False)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'SUTM1', help='Sparse upper triangular matrix file, group 1.')
    parser.add_argument(
        'SUTM2', help='Sparse upper triangular matrix file, group 2.')
    parser.add_argument(
        '--seed', default=None, type=int,
        help='Seed for random swap (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--out1', required=True,
        help='Path prefix of shadow matrices (group 1) in SUTM format.')
    requiredNamed.add_argument(
        '--out2', required=True,
        help='Path prefix of shadow matrices (group 2) in SUTM format.')
    parser.set_defaults(function=homer2sutm)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
