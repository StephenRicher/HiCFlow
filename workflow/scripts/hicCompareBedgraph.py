#!/usr/bin/env python3

""" Compute sum logFC for all, up and down interactions """

import sys
import argparse
import pandas as pd
from scipy.stats import zscore
from utilities import setDefaults, createMainParent, readHomer


__version__ = '1.0.0'


def hicCompareBedgraph(
        file: str, allOut: str, upOut: str, downOut: str,
        binSize: int, maxDistance: float):

    positions, mat = readHomer(file, binSize)
    mat['seperation'] = abs(mat['end'] - mat['start'])
    # Remove paired interactions above max distance
    if maxDistance:
        mat = mat.loc[mat['seperation'] < maxDistance]
    mat['upFC'] = mat['score'] > 0
    mat['score'] = abs(mat['score'])
    mat = mat.set_index('region').loc[:,['score', 'upFC']]

    for direction in ['up', 'down', 'all']:
        if direction == 'up':
            out = upOut
            subset = mat.loc[mat['upFC'] == True, 'score'].groupby('region').sum()
        elif direction == 'down':
            out = downOut
            subset = mat.loc[mat['upFC'] == False, 'score'].groupby('region').sum()
        else:
            out = allOut
            subset = mat.loc[:, 'score'].groupby('region').sum()

        bed = pd.merge(
            positions, subset, how='left',
            left_index=True, right_index=True).fillna(0)
        bed['zscore'] = zscore(bed['score'])
        bed.to_csv(
            out,  columns=[0, 1, 2, 'zscore'],
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
