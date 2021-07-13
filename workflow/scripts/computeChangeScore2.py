#!/usr/bin/env python3

""" Compute permutation score of between sum(logFC) at bin position """

import os
import sys
import argparse
import numpy as np
import pandas as pd
from typing import List
from matplotlib import cm
from functools import reduce
from matplotlib.colors import to_hex
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def computeChangeScore(matrix: str, shadowMatrices: List, rawOut: str, binSize: int, chrom: str):

    allRegions = []
    nShadow = len(shadowMatrices)
    for i, mat in enumerate([matrix] + shadowMatrices):
        mat = readSUTM(mat, lower=True)
        mat[i] = abs(mat['score'])
        summed = mat.groupby('start1')[[i]].sum().reset_index()
        allRegions.append(summed.set_index('start1'))
    allRegions = reduce(
        lambda l, r: pd.merge(
            l, r, left_index=True, right_index=True, how='left'), allRegions)

    allRegions = (allRegions
        .fillna(0)
        .melt(var_name='shadow', value_name='abs(score)', ignore_index=False)
        .reset_index())
    allRegions['real'] = allRegions['shadow'] == 0

    permutation = (allRegions.groupby('start1')
        .apply(permutationTest, nShadow=nShadow)
        .reset_index().rename({0: 'p'}, axis=1))

    permutation['chrom'] = chrom
    permutation['name'] = '.'
    permutation['end'] = permutation['start1'] + binSize
    columns = ['chrom', 'start1', 'end', 'name', 'p']
    permutation[columns].to_csv(sys.stdout, index=False, header=False, sep='\t')

    if rawOut is not None:
        allRegions['chrom'] = chrom
        allRegions.to_csv(rawOut, header=True, sep='\t')


def readSUTM(sutm, lower=False):
    sutm = pd.read_csv(sutm, names=['start1', 'start2', 'score'], sep='\s+')
    if lower:
        sltm = sutm.loc[sutm['start1'] != sutm['start2']].rename(
            {'start1': 'start2', 'start2': 'start1'}, axis=1)
        sutm  = pd.concat([sutm, sltm])
    return sutm


def permutationTest(x, nShadow):
    """ Compare shadow abs(score) against normal matrix """
    shadowScores = x.loc[~x['real'], 'abs(score)']
    normalScore = float(x.loc[x['real'], 'abs(score)'])
    totalAbove = (shadowScores >= normalScore).sum()
    return (totalAbove / nShadow) + (1 / nShadow)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=computeChangeScore)
    parser.add_argument(
        'matrix', help='HiC matrix in SUTM format.')
    parser.add_argument(
        'shadowMatrices', nargs='*',
        help='Shadow HiC matrices inSUTM format.')
    parser.add_argument(
        '--chrom', help='Reference of correponding matrices.')
    parser.add_argument(
        '--binSize', type=int, help='Bin size of corresponding matrices')
    parser.add_argument(
        '--rawOut', help='Path for raw, unsummarised data.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
