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


def readHomer(matrix, binSize=None, sparse=True):
    """ Read Homer matrix format and convert to long format """

    # Read Homer as pandas
    mat = pd.read_csv(matrix, skiprows=1, header=None, sep='\t').drop(0, axis=1)
    # Split chromosome and start into 2 column DF
    regions = mat[1].str.split('-', expand=True)
    startPos = regions[1].astype(int)
    # Drop region column
    mat = mat.drop(1, axis=1)
    # Set columns as bin start positions
    mat.columns = list(startPos)
    # Set rownames (index) as bin start positions
    mat['start'] = list(startPos)
    mat = mat.set_index('start')
    # Convert to long format
    mat = mat.stack().rename_axis(['start', 'start2']).rename('score').reset_index()
    # Remove 0 score rows if sparse set
    if sparse:
        mat = mat.loc[mat['score'] != 0]
    # Remove symetric interactions
    mat = mat.loc[mat['start2'] - mat['start'] >= 0]
    # Set chrom from first value since HOMER must be cis-only matrix
    mat.attrs['chrom'] = regions[0][0]
    mat.attrs['binSize'] = binSize

    return mat
