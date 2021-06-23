#!/usr/bin/env python3

""" Compute change score and estimate directional preference """

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


def computeASAP(matrices: List, outData: str, colourmap: str, fdr: float):

    scalingFactor = computeScalingFactor(matrices)
    allRegions = []
    allDirection = []
    for matrix in matrices:
        mat = readHomer(matrix, diagonal=True, sparse=False)
        chrom = mat.attrs['chrom']
        binSize = mat.attrs['binSize']
        binDirection = directionPreference(mat, chrom, binSize)
        allDirection.append(binDirection)
        mat['seperation'] = (mat['start'] - mat['start2']).abs()
        mat = pd.merge(mat, scalingFactor, left_on='seperation', right_on='seperation')
        mat['score'] = mat['score'] * mat['scale']
        mat['abs(score)'] = abs(mat['score'])
        summed = mat.groupby('start')[['abs(score)']].sum().reset_index()
        summed['chrom'] = chrom
        allRegions.append(summed)

    # Combine and shuffle rows to avoid bias in ranking
    allRegions = pd.concat(allRegions).set_index(['chrom', 'start']).sample(frac=1)
    allRegions['quantScore'] = pd.qcut(
        allRegions['abs(score)'].rank(method='first'), 20, labels=np.linspace(0, 1, 20))
    allDirection = pd.concat(allDirection).set_index(['chrom', 'start'])

    allRegions = allRegions.merge(
        allDirection, left_index=True, right_index=True).reset_index()
    allRegions['end'] = allRegions['start'] + binSize
    validP = allRegions['p'].notna()
    allRegions['p(adj)'] = np.nan
    allRegions.loc[validP, 'p(adj)'] = fdrcorrection(allRegions.loc[validP, 'p'])[1]

    #allRegions['p(adj)'] = fdrcorrection(allRegions['p'])[1]
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


def computeScalingFactor(matrices: List):
    """ Read all matrices to obtain for range of bin seperations.
        Apply a power law scaling factor between 0 and 1 to downweight
        more distant interactions """

    allBinSizes = set()
    allSeperations = set()
    for matrix in matrices:
        with open(matrix) as fh:
            header = fh.readline().strip().split('\t')[2:]
            header = pd.Series([int(pos.split('-')[1]) for pos in header])
            try:
                binSize = int(header.diff().dropna().unique())
            except TypeError:
                sys.exit(f'Bin sizes in {matrix} are not all equal.')
            allBinSizes.add(binSize)
            binRange = header.iloc[-1] - header.iloc[0]
            allSeperations.update(list(range(0, binRange + 1, binSize)))
    assert len(allBinSizes) == 1, 'Binsizes must be equal across matrices'
    allSeperations = list(sorted(allSeperations))
    scalingFactor = np.power(np.linspace(1.0, 0.0, num=len(allSeperations)), 1)
    scalingFactor = pd.DataFrame(
        {'seperation': allSeperations, 'scale': scalingFactor})

    return scalingFactor


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
    parser.set_defaults(function=computeASAP)
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


    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
