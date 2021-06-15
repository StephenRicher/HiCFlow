#!/usr/bin/env python3

""" Compute absolute logFC and discretise """

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from typing import List
from matplotlib import cm
from matplotlib.colors import to_hex
from collections import defaultdict
from statsmodels.stats.multitest import fdrcorrection
from utilities import setDefaults, createMainParent, readHomer


__version__ = '1.0.0'


def shadowCompareHiC(
        matrices: List, outData: str, colourmap: str,
        fdr: float, maxDistance: int):

    allRegions = []
    allDirection = []
    for matrix in matrices:
        mat = readHomer(matrix, diagonal=True, sparse=True)
        chrom = mat.attrs['chrom']
        binSize = mat.attrs['binSize']

        binDirection = directionPreference(mat, chrom, binSize)
        allDirection.append(binDirection)

        # Find and drop bins with a few non-missing values
        binZeros = mat.groupby('start')['score'].count()
        dropBins = binZeros[binZeros < 10].index
        dropRows = mat['start'].isin(dropBins)
        mat = mat.loc[~dropRows]

        # Remove interactions above a certain distance
        mat['seperation'] = (mat['start'] - mat['start2']).abs()
        mat = mat.loc[abs(mat['seperation']) < maxDistance]

        mat['abs(score)'] = abs(mat['score'])
        summed = mat.groupby('start')[['abs(score)']].mean().reset_index()

        summed['chrom'] = chrom
        allRegions.append(summed)

    allRegions = pd.concat(allRegions).set_index(['chrom', 'start'])
    allRegions['quantScore'] = pd.qcut(
        allRegions['abs(score)'], 20, labels=np.linspace(0, 1, 20))
    allDirection = pd.concat(allDirection).set_index(['chrom', 'start'])

    allRegions = allRegions.merge(
        allDirection, left_index=True, right_index=True).reset_index()
    allRegions['end'] = allRegions['start'] + binSize
    allRegions['p(adj)'] = fdrcorrection(allRegions['p'])[1]
    allRegions['colour'] = allRegions.apply(getColour, args=(fdr, colourmap), axis=1)

    if outData is not None:
        allRegions.to_csv(
            outData, index=False, header=True, sep='\t')

    # Write direction score as a BED file
    allRegions['name'] = '.'
    allRegions['strand'] = '.'
    allRegions['thickStart'] = allRegions['start']
    allRegions['thickEnd'] = allRegions['end']
    columns = ([
        'chrom', 'start', 'end', 'name', 'abs(score)', 'strand',
        'thickStart', 'thickEnd', 'colour'])
    allRegions[columns].to_csv(
        sys.stdout, index=False, header=False, sep='\t')



def getColour(x, fdr, colourmap):
    if x['p(adj)'] <= fdr:
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
        _, p = stats.wilcoxon(df['score'], alternative='two-sided')
        diffChange['p'].append(p)
        direction = 1 if df['score'].median() > 0 else -1
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
    parser.set_defaults(function=shadowCompareHiC)
    parser.add_argument(
        'matrices', nargs='*', help='HiC matrix in homer format.')
    parser.add_argument(
        '--outData', default=None,
        help='File to write all relevant difference score information '
             '(default: %(default)s)')
    parser.add_argument(
        '--fdr', type=float, default=0.05,
        help='False discovery rate threshold for directional '
             'bias sigificance (default: %(default)s)')
    parser.add_argument(
        '--colourmap', default='bwr',
        help='Colour map for directional change (default: %(default)s)')
    parser.add_argument(
        '--maxDistance', type=int, default=1000000,
        help='Remove interactions greater than this distance '
             '(default: %(default)s)')


    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
