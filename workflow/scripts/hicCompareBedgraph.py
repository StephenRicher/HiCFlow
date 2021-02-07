#!/usr/bin/env python3

""" Compute sum logFC for all, up and down interactions """

import sys
import argparse
import pandas as pd
from scipy import stats
from utilities import setDefaults, createMainParent, readHomer


__version__ = '1.0.0'


def hicCompareBedgraph(
        file: str, allOut: str, upOut: str = None,
        downOut: str = None, minDistance: float = None,
        maxDistance: float = None, Z: bool = False):

    mat = readHomer(file, diagonal=False, sparse=True)
    # Retrieve all matrix start positions
    allStart = pd.Series(mat['start'].unique(), name='start')
    # Remove paired interactions above max distance
    if maxDistance:
        inRange = abs(mat['start2'] - mat['start']) < abs(maxDistance)
        mat = mat.loc[inRange]
    if minDistance:
        inRange = abs(mat['start2'] - mat['start']) > abs(minDistance)
        mat = mat.loc[inRange]
    mat['upFC'] = mat['score'] > 0
    mat['score'] = abs(mat['score'])

    mat = mat.set_index('start').loc[:,['score', 'upFC']]

    for direction in ['up', 'down', 'all']:
        if direction == 'up':
            if upOut is None:
                continue
            out = upOut
            subset = mat.loc[mat['upFC'] == True, 'score'].groupby('start').sum()
        elif direction == 'down':
            if downOut is None:
                continue
            out = downOut
            subset = mat.loc[mat['upFC'] == False, 'score'].groupby('start').sum()
        else:
            out = allOut
            subset = mat.loc[:, 'score'].groupby('start').sum()

        if Z:
            subset = stats.zscore(subset)
        score = pd.Series(subset, index=subset.index, name='score')
        bed = pd.merge(
            allStart, zscore, how='left', left_on='start', right_index=True).fillna(0)
        bed['chrom'] = mat.attrs['chrom']
        bed['end'] = bed['start'] + mat.attrs['binSize']

        bed.to_csv(
            out,  columns=['chrom', 'start', 'end', 'score'],
            sep='\t', index=False, header=False)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=hicCompareBedgraph)
    parser.add_argument('file', help='HiC matrix in homer format.')
    parser.add_argument(
        '--minDistance', type=float,
        help='Only consider interactions above this distance (default: %(default)s)')
    parser.add_argument(
        '--maxDistance', type=float,
        help='Only consider interactions less this distance.')
    parser.add_argument(
        '--upOut',
        help='Output file for bedgraph of sum UP interactions.')
    parser.add_argument(
        '--downOut',
        help='Output file for bedgraph of sum down interactions.')
    parser.add_argument(
        '--Z', action='store_true',
        help='Z score transform scores (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--allOut', required=True,
        help='Output file for bedgraph of sum all interactions.')


    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
