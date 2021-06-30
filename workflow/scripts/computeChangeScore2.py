#!/usr/bin/env python3

""" Compute change score and estimate directional preference """

import os
import sys
import argparse
import pandas as pd
from typing import List
from matplotlib import cm
from matplotlib.colors import to_hex
from utilities import setDefaults, createMainParent, readHomer


__version__ = '1.0.0'


def computeChangeScore(matrix: str, shadowMatrices: List):

    allRegions = []
    for i, mat in enumerate([matrix] + shadowMatrices):
        mat = readHomer(mat, diagonal=True, sparse=False)
        chrom = mat.attrs['chrom']
        binSize = mat.attrs['binSize']
        mat['seperation'] = (mat['start'] - mat['start2']).abs()
        mat['abs(score)'] = abs(mat['score'])
        summed = mat.groupby('start')[['abs(score)']].sum().reset_index()
        summed['chrom'] = chrom
        if i == 0:
            summed['shadow'] = False
        else:
            summed['shadow'] = True
        allRegions.append(summed)

    # Combine and shuffle rows to avoid bias in ranking
    allRegions = pd.concat(allRegions)
    permutation = (allRegions.groupby(['chrom', 'start']).apply(permutationTest)
        .reset_index().rename({0: 'p'}, axis=1))

    permutation['score'] = 1 - permutation['p']
    permutation['colour'] = permutation.apply(getColour, axis=1)
    permutation['name'] = '.'
    permutation['strand'] = '.'
    permutation['end'] = permutation['start'] + binSize
    permutation['thickStart'] = permutation['start']
    permutation['thickEnd'] = permutation['end']
    columns = ([
        'chrom', 'start', 'end', 'name', 'score', 'strand',
        'thickStart', 'thickEnd', 'colour'])
    permutation[columns].to_csv(sys.stdout, index=False, header=False, sep='\t')


def permutationTest(x):
    """ Compare shadow abs(score) against normal matrix """
    shadowScores = x.loc[x['shadow'], 'abs(score)']
    normalScore = float(x.loc[~x['shadow'], 'abs(score)'])
    totalAbove = (shadowScores > normalScore).sum()
    return totalAbove / len(shadowScores)


def getColour(x):
    colour = to_hex(cm.get_cmap('binary', 100)(x['score']))[1:]
    colour = f'{int(colour[:2], 16)},{int(colour[2:4], 16)},{int(colour[4:], 16)}'
    return colour


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=computeChangeScore)
    parser.add_argument('matrix', help='HiC matrix in homer format.')
    parser.add_argument(
        'shadowMatrices', nargs='*',
        help='Shadow HiC matrices in homer format.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
