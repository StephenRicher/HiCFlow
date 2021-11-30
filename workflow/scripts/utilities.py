#!/usr/bin/env python3

import os
import sys
import logging
import argparse
import itertools
import contextlib
import numpy as np
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


class StoreDict(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        kv={}
        if not isinstance(values, (list,)):
            values=(values,)
        for value in values:
            if value.count('=') != 1:
                logging.error(
                    f'Key value data in {value} cannot contain exactly one "=".')
                raise ValueError
            n, v = value.split('=')
            kv[n]=v
        setattr(namespace, self.dest, kv)


def readHomer(matrix, upperOnly=False, sparse=True, diagonal=False, distanceNorm=False, absolute=False):
    """ Read Homer matrix format and convert to long format """

    # Read Homer as pandas
    mat = pd.read_csv(matrix, skiprows=1, header=None, sep='\t').drop(0, axis=1)
    # Split chromosome and start into 2 column DF
    regions = mat[1].str.split('-', expand=True)
    startPos = regions[1].astype(int)
    # Ensure homer is constant binSize matrix
    binSize = startPos.diff(1).dropna().unique()
    assert len(binSize) == 1
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
    if absolute:
        mat['score'] = mat['score'].abs()
    if distanceNorm:
        mat['seperation'] = abs(mat['start2'] - mat['start'])
        sumDistance = mat.groupby('seperation').score.agg(['sum', 'count'])
        sumDistance['expected'] = sumDistance['sum'] / sumDistance['count']
        mat = pd.merge(mat, sumDistance, on="seperation")
        mat['score'] = mat['score'] / mat['expected']
    if upperOnly:
        mat = mat.loc[mat['start2'] > mat['start']]
    if diagonal is False:
        mat = mat.loc[mat['start2'] != mat['start']]
    # Set chrom from first value since HOMER must be cis-only matrix
    mat.attrs['chrom'] = regions[0][0]
    mat.attrs['binSize'] = int(binSize)

    return mat[['start', 'start2', 'score']]


def readChromSizes(file):
    """ Read chromosome sizes to dict """
    chromSizes = {}
    with open(file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chrom, size = line.split()
            chromSizes[str(chrom)] = int(size)
    return chromSizes


def getGroup(path):
    """ Retrieve group name from file path """
    basename = os.path.basename(path)
    return basename.split('-')[0]


def getBin(pos: int, binSize: int, start: int):
    matrixBin = (pos - start) // binSize
    return start + (matrixBin * binSize)


def readSUTM(sutm, lower=False, bothDiagonal=False):
    sutm = pd.read_csv(sutm, names=['start1', 'start2', 'adjIF'], sep=' ')
    sutm['orientation'] = 1
    if lower:
        # Remove the diagonal for the opposite orientation
        if not bothDiagonal:
            sltm = sutm.loc[sutm['start1'] != sutm['start2']]
        else:
            sltm = sutm.copy()
        sltm = sltm.rename({'start1': 'start2', 'start2': 'start1'}, axis=1)
        sltm['orientation'] = -1
        sutm  = pd.concat([sutm, sltm])
    return sutm.set_index(['start1', 'start2', 'orientation'])
