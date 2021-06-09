#!/usr/bin/env python3

""" Convert bedgraph of subtraction scores to bed of regions with > 0 Z """

import sys
import logging
import argparse
import itertools
import pandas as pd
from scipy.stats import zscore
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def subtraction2Bed(bedgraph: str, bedgraph2: str):
    dtypes = {'chrom': str, 'start': int, 'end': int, 'score': float}
    bedgraph = pd.read_csv(
        bedgraph, names=dtypes.keys(), dtype=dtypes, sep='\t')
    bedgraph['name'] = '.'
    if bedgraph2:
        bedgraph2 = pd.read_csv(
            bedgraph2, names=dtypes.keys(), dtype=dtypes, sep='\t')
        # Ensure bedgraph coordinates match
        headers = ['chrom', 'start', 'end']
        assert bedgraph[headers].equals(bedgraph2[headers])
        bedgraph['score2'] = bedgraph2['score']
        bedgraph['score1higher'] = bedgraph.apply(
            lambda x: 1 if x['score'] > x['score2'] else -1, axis=1)
        bedgraph = bedgraph.loc[bedgraph[['score','score2']].max(axis=1) > 0]
    else:
        bedgraph = bedgraph.loc[bedgraph['score'] > 0]
        bedgraph['score1higher'] = 1
    bedgraph[['chrom', 'start', 'end', 'name', 'score1higher']].to_csv(
        sys.stdout, header=False, index=False, sep='\t')

def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'bedgraph',
        help='Bedgraph subtraction score.')
    parser.add_argument(
        'bedgraph2',  nargs='?',
        help='Bedgraph subtraction score to compare against.')
    parser.set_defaults(function=subtraction2Bed)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
