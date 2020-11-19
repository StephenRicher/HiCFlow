#!/usr/bin/env python3

""" Compute sum logFC for up and down interactions """

import sys
import argparse
import pandas as pd
from utilities import setDefaults

__version__ = '1.0.0'


def hicCompareBedgraph(file: str, upOut: str, downOut: str, binSize: int):

    positions, mat = readHomer(file, binSize)
    mat['upFC'] = mat['score'] > 0
    mat['score'] = abs(mat['score'])
    mat = mat.groupby(['bin', 'upFC']).sum().reset_index('upFC')

    for upFC in [True, False]:
        out = upOut if upFC else downOut
        bed = pd.merge(
            positions, mat.loc[mat['upFC'] == upFC],
            how='left', left_index=True, right_index=True)
        bed['zscore'] = (bed['score'] - bed['score'].mean()) / bed['score'].std(ddof=0)
        bed['zscore'] = bed['zscore'].fillna(0)
        bed.to_csv(
            out,  columns=[0, 1, 2, 'zscore'],
            sep='\t', index=False, header=False)


def readHomer(matrix, binSize):
    """ Read Homer matrix format and convert to long format """
    mat = pd.read_csv(matrix, skiprows=1, header=None, sep='\t').drop(0, axis=1)
    # Extract chromsome and start coordinates
    positions = mat[1].str.split('-', expand=True)
    # Add end positions
    positions[2] = positions[1].astype(int) + binSize
    positions['pos'] = mat[1]
    positions =  positions.set_index('pos')
    # Convert to long and drop 'otherBin' column
    mat = mat.melt(id_vars=1).drop('variable', axis=1)
    mat.columns = ['bin', 'score']

    return positions, mat


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('file', help='HiC matrix in homer format.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--binSize', required=True,
        type=int, help='Bin size for matrix.')
    requiredNamed.add_argument(
        '--upOut', required=True,
        help='Output file for bedgraph of sum UP interactions.')
    requiredNamed.add_argument(
        '--downOut', required=True,
        help='Output file for bedgraph of sum down interactions.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(hicCompareBedgraph(**vars(args)))
