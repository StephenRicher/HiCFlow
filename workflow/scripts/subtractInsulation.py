#!/usr/bin/env python3

""" Compare insulation score between a pair of insulation scores """

import sys
import logging
import argparse
import itertools
import pandas as pd
from scipy.stats import zscore
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def subtractInsulation(TADinsulation1: str, TADinsulation2: str):
    TADinsulation1 = pd.read_csv(
        TADinsulation1, index_col=[0, 1, 2], header=None, comment='#', sep='\t')
    TADinsulation2 = pd.read_csv(
        TADinsulation2, index_col=[0, 1, 2], header=None, comment='#', sep='\t')

    # Perform difference on each TAD insulation and compute absolute difference
    TADdifference = (TADinsulation1.diff(1) - TADinsulation2.diff(1)).abs().fillna(0).apply(zscore)
    TADdifference.to_csv(sys.stdout, header=False, sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'TADinsulation1', help='TAD insulation score as bedgraph matrix.')
    parser.add_argument(
        'TADinsulation2', help='TAD insulation score as bedgraph matrix.')
    parser.set_defaults(function=subtractInsulation)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
