#!/usr/bin/env python3

""" Compute sum logFC for all, up and down interactions """

import sys
import argparse
import pandas as pd
from scipy import stats
from utilities import setDefaults, createMainParent, readHomer


__version__ = '1.0.0'


def hicCompareBedgraph(
        file: str, allOut: str, upOut: str, downOut: str,
        binSize: int, maxDistance: float):

    mat = readHomer(file, binSize, sparse=True)
    # Retrieve all matrix start positions
    allStart = pd.Series(mat['start'].unique(), name='start')
    # Remove paired interactions above max distance
    if maxDistance:
        inRange = abs(mat['start2'] - mat['start']) < maxDistance
        mat = mat.loc[inRange]
    mat['upFC'] = mat['score'] > 0
    mat['score'] = abs(mat['score'])

    mat = mat.set_index('start').loc[:,['score', 'upFC']]

    for direction in ['up', 'down', 'all']:
        if direction == 'up':
            out = upOut
            subset = mat.loc[mat['upFC'] == True, 'score'].groupby('start').sum()
        elif direction == 'down':
            out = downOut
            subset = mat.loc[mat['upFC'] == False, 'score'].groupby('start').sum()
        else:
            out = allOut
            subset = mat.loc[:, 'score'].groupby('start').sum()

        # Perform Z-score normalisation
        subset.update(pd.Series(stats.zscore(subset)))

        bed = pd.merge(
            allStart, subset, how='left', left_on='start', right_index=True).fillna(0)
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
        '--maxDistance', type=float,
        help='Only consider interactions within this distance.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--binSize', required=True,
        type=int, help='Bin size for matrix.')
    requiredNamed.add_argument(
        '--allOut', required=True,
        help='Output file for bedgraph of sum all interactions.')
    requiredNamed.add_argument(
        '--upOut', required=True,
        help='Output file for bedgraph of sum UP interactions.')
    requiredNamed.add_argument(
        '--downOut', required=True,
        help='Output file for bedgraph of sum down interactions.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
