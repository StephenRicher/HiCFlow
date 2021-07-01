#!/usr/bin/env python3

""" Calculate adjusted absolute logFC from adjIF1 and adjIF2 SUTM matrices
    computed by HiCcompare. Also calculate per-bin merged counts for procesing
    by shadowSUTM.py. """

import os
import sys
import argparse
import numpy as np
import pandas as pd
from typing import List
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def computeAdjM(adjIF1: str, adjIF2: str, adjM_out: str, merge_out: str):

    adjIF1 = readSUTM(adjIF1, lower=True, integer=True)
    adjIF2 = readSUTM(adjIF2, lower=True, integer=True)

    # Compute absolute logFC
    adjM = pd.merge(adjIF1, adjIF2, left_index=True, right_index=True)
    adjM[0] = np.log2(adjM['adjIF_x'] / adjM['adjIF_y'])#.abs()

    # Remove bins with partial counts
    adjM = adjM.loc[~adjM[0].isin([np.nan, np.inf, -np.inf]),]

    # Retrieve ratio of IF between samples to set sampling probability
    totalIF = adjIF1['adjIF'].sum() + adjIF2['adjIF'].sum()
    IF1ratio = adjIF1['adjIF'].sum() / totalIF

    # Sum values across bins and save to pickle
    adjM = adjM.groupby('start1')[0].sum().to_frame().to_pickle(adjM_out)

    # Combine bin-pair counts of both groups
    mergedSUTM = pd.concat([adjIF1, adjIF2]).groupby(['start1', 'start2']).sum()
    mergedSUTM.attrs['IF1ratio'] = IF1ratio
    mergedSUTM.to_pickle(merge_out)


def readSUTM(sutm, lower=False, integer=False):
    sutm = pd.read_csv(sutm, names=['start1', 'start2', 'adjIF'], sep=' ')
    if lower:
        sltm = sutm.loc[sutm['start1'] != sutm['start2']].rename(
            {'start1': 'start2', 'start2': 'start1'}, axis=1)
        sutm  = pd.concat([sutm, sltm])
    if integer:
        sutm['adjIF'] = sutm['adjIF'].astype(int)
    return sutm.set_index(['start1', 'start2'])


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
    parser.set_defaults(function=computeAdjM)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
