#!/usr/bin/env python3

""" Compute change score and estimate directional preference """

import os
import sys
import argparse
import numpy as np
import pandas as pd
from typing import List
from scipy import stats
from matplotlib import cm
from collections import defaultdict
from matplotlib.colors import to_hex
from hicmatrix import HiCMatrix as hm
from sklearn.preprocessing import KBinsDiscretizer
from statsmodels.stats.multitest import fdrcorrection
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def computeChangeScore(matrices: List, colourmap: str, alpha: float,
                       maxDistance: int, nBins: int):

    allRegions = []
    allDirection = []
    est = KBinsDiscretizer(n_bins=nBins, encode='ordinal', strategy='kmeans')
    for matrix in matrices:
        hic = hm.hiCMatrix(matrix)
        binSize = hic.getBinSize()
        chrom = hic.getChrNames()[0]

        nonzeroIdx = hic.matrix.nonzero()
        nonzeroValues = hic.matrix[nonzeroIdx].tolist()[0]
        mat = pd.DataFrame({
            'start' : nonzeroIdx[0],
            'start2': nonzeroIdx[1],
            'value' : nonzeroValues
        })
        # Associtate bin ID with genomic position
        mat['start'] = mat['start'].apply(lambda x: hic.getBinPos(x)[1])
        mat['start2'] = mat['start2'].apply(lambda x: hic.getBinPos(x)[1])
        mat.columns = ['start', 'start2', 'score']
        mat['seperation'] = (mat['start'] - mat['start2']).abs()
        mat['abs(score)'] = abs(mat['score'])
        mat = mat.loc[mat['seperation'] <= maxDistance]

        binDirection = directionPreference(mat, chrom, binSize)

        # FDR correction by chromosome
        validP = binDirection['p'].notna()
        binDirection['p(adj)'] = np.nan
        binDirection.loc[validP, 'p(adj)'] = fdrcorrection(
            binDirection.loc[validP, 'p'])[1]
        allDirection.append(binDirection)

        summed = mat.groupby('start')[['abs(score)']].sum().reset_index()
        summed['chrom'] = chrom
        summed['quantScore'] = est.fit_transform(summed['abs(score)'].to_numpy().reshape(-1,1)) / nBins
        allRegions.append(summed)

    allRegions = pd.concat(allRegions).set_index(['chrom', 'start'])
    allDirection = pd.concat(allDirection).set_index(['chrom', 'start'])

    allRegions = allRegions.merge(
        allDirection, left_index=True, right_index=True).reset_index()
    allRegions['end'] = allRegions['start'] + binSize

    allRegions['colour'] = allRegions.apply(
        getColour, args=(alpha, colourmap), axis=1)

    allRegions['name'] = '.'
    allRegions['strand'] = '.'
    allRegions['thickStart'] = allRegions['start']
    allRegions['thickEnd'] = allRegions['end']
    columns = ([
        'chrom', 'start', 'end', 'name', 'abs(score)', 'strand',
        'thickStart', 'thickEnd', 'colour'])
    allRegions[columns].to_csv(sys.stdout, index=False, header=False, sep='\t')


def getColour(x, alpha, colourmap, p='p(adj)'):
    if x[p] <= alpha:
        if x['direction'] == 1:
            i = (x['quantScore'] * 0.5)  + 0.5
        else:
            i = (1 - x['quantScore']) * 0.5
        colour = to_hex(cm.get_cmap(colourmap, 40)(i))[1:]
    else:
        i = x['quantScore']
        colour = to_hex(cm.get_cmap('binary', 20)(i))[1:]
    colour = f'{int(colour[:2], 16)},{int(colour[2:4], 16)},{int(colour[4:], 16)}'
    return colour


def directionPreference(matrix, chrom: str, binSize: int):

    diffChange = defaultdict(list)
    for start, df in matrix.groupby('start'):
        diffChange['start'].append(start)
        try:
            _, p = stats.wilcoxon(df['score'], alternative='two-sided')
        except ValueError:
            p = np.nan

        diffChange['p'].append(p)
        nonZeroScores = df.loc[df['score'] != 0, 'score']
        if not nonZeroScores.empty:
            direction = 1 if nonZeroScores.median() > 0 else -1
        else:
            direction = np.nan
        diffChange['direction'].append(direction)
    diffChange = pd.DataFrame(diffChange)
    diffChange['chrom'] = chrom
    diffChange['end'] = diffChange['start'] + binSize

    return diffChange


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=computeChangeScore)
    parser.add_argument(
        'matrices', nargs='+', help='HiC matrix in homer format.')
    parser.add_argument(
        '--alpha', type=float, default=0.05,
        help='False discovery rate threshold for directional '
             'bias sigificance (default: %(default)s)')
    parser.add_argument(
        '--colourmap', default='bwr',
        help='Colour map for directional change (default: %(default)s)')
    parser.add_argument(
        '--maxDistance', type=int, default=1e6,
        help='Maximum interaction distance (bp) to compute '
             'change score (default: %(default)s)')
    parser.add_argument(
        '--nBins', type=int, default=20,
        help='Number of clusters to group change '
             'score (default: %(default)s)')


    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
