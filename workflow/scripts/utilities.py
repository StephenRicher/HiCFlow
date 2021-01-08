#!/usr/bin/env python3

import sys
import logging
import argparse
import pandas as pd
from collections import defaultdict


def setDefaults(parser):
    """ Set logging config and return args with associated function """

    args = parser.parse_args()
    logFormat='%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
    try:
        logging.basicConfig(level=args.verbose, format=logFormat)
        del args.verbose
    except AttributeError:
        logging.basicConfig(format=logFormat)
        pass
    try:
        function = args.function
        del args.function
    except AttributeError:
        parser.print_help()
        sys.exit()

    return args, function


def createMainParent(verbose=True, version=None):
    """ Create parser of verbose/version to be added to parser/subparsers. """
    parent = argparse.ArgumentParser(add_help=False)
    if version:
        parent.add_argument('--version', action='version',
            version=f'%(prog)s {version}')
    if verbose:
        parent.add_argument(
            '--verbose', action='store_const', const=logging.DEBUG,
            default=logging.INFO, help='verbose logging for debugging')
    return parent


def readHomer(matrix, binSize):
    """ Read Homer matrix format and convert to long format """

    # Read Homer as pandas
    mat = pd.read_csv(matrix, skiprows=1, header=None, sep='\t').drop(0, axis=1)

    # Split chromosome and start into 2 column DF
    positions = mat[1].str.split('-', expand=True)
    positions[1] = positions[1].astype(int)
    # Add end positions
    positions[2] = positions[1] + binSize
    # Add full region position and set as index
    positions['pos'] = mat[1]
    positions =  positions.set_index('pos')

    # Drop region column
    mat = mat.drop(1, axis=1)
    # Set index and columns as start positions
    mat.columns = list(positions[1])

    mat['start'] = mat.columns
    mat = mat.melt(id_vars=['start'], var_name='start2', value_name='score')
    mat['end'] = mat['start'] + binSize
    mat['end2'] = mat['start2'] + binSize

    # Extract chrom - assume matrix must be cis
    chrom = f'{positions[0][0]}'
    mat.attrs['chrom'] = chrom
    mat['region'] = f'{chrom}-' + mat['start'].astype(str)

    return positions, mat
