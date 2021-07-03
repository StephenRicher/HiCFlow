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
from statsmodels.stats.multitest import fdrcorrection
from utilities import setDefaults, createMainParent
try:
    from pandarallel import pandarallel
except ModuleNotFoundError:
    pass


__version__ = '1.0.0'


def computeAdjM(
        adjIF1: str, adjIF2: str, fdr: float, nShadow: int,
        rawOut: str, binSize: int, chrom: str, threads: int):

    adjIF1 = readSUTM(adjIF1, lower=True)
    adjIF2 = readSUTM(adjIF2, lower=True)

    # Compute absolute logFC
    adjM = pd.merge(
        adjIF1, adjIF2, how='inner',
        left_index=True, right_index=True)
    adjM = adjM.loc[adjM['adjIF_x'] != adjM['adjIF_y']]

    adjM['diff'] = (adjM['adjIF_x'] - adjM['adjIF_y']).apply(getDirection)

    if (threads > 1) and ('pandarallel' in sys.modules):
        pandarallel.initialize(nb_workers=5, verbose=0)
        permuted = adjM.groupby('start1').parallel_apply(permuteTest, nShadow)
    else:
        permuted = adjM.groupby('start1').apply(permuteTest, nShadow)
    permuted = permuted.reset_index().rename({0: 'p'}, axis=1)

    permuted['fdr'] =  fdrcorrection(permuted['p'])[1]
    permuted['chrom'] = chrom
    permuted['end'] = permuted['start1'] + binSize
    columns = ['chrom', 'start1', 'end', 'p', 'fdr']
    if rawOut:
        permuted.to_csv(rawOut, index=False, header=False, sep='\t')
    permuted = permuted.loc[permuted['fdr'] < fdr]
    permuted[columns].to_csv(sys.stdout, index=False, header=False, sep='\t')


def permuteTest(x, nShadow):
    length = rle(np.array([x['diff'].values]))
    boolRand = np.random.random((nShadow, len(x))) > 0.5
    randRLE = rle(boolRand)
    totalAbove = (randRLE <= length).sum()
    p = (totalAbove / nShadow)
    if p == 0:
        p += (1 / nShadow)
    return p


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
        y = ia[...,1:] != ia[..., :-1]               # pairwise unequal (string safe)
        a = np.append(y, np.ones((len(y), 1)), 1)
        x = np.where(a)[1]
        indexes = np.where(x == y.shape[1])[0] + 1
        splits = np.split(x, indexes)[:-1]
        return np.array([len(np.diff(x)) + 1 for x in splits])


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
        '--fdr', default=0.1, type=float,
        help='FDR significance threshold (default: %(default)s)')
    parser.add_argument(
        '--rawOut', help='Output file for unfiltered results')
    parser.add_argument(
        '--nShadow', default=1000, type=int,
        help='Number of random permutations (default: %(default)s)')
    parser.add_argument(
        '--threads', default=1, type=int,
        help='Threads for parallel processing (default: %(default)s)')
    parser.add_argument(
        '--chrom', help='Reference of correponding matrices.')
    parser.add_argument(
        '--binSize', type=int, help='Bin size of corresponding matrices')
    parser.set_defaults(function=computeAdjM)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
