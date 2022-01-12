#!/usr/bin/env python3

""" Correlate Cscore and write absolute difference to BED """


import sys
import argparse
import pandas as pd
from scipy.stats import pearsonr
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def computeSwitchScore(cscore1: str, cscore2: str) -> None:

    names = {'chrom': str, 'start': int, 'end': int, 'cscore': float}
    cscore1 = pd.read_csv(
        cscore1, usecols=[0,1,2,4], names=names.keys(),
        dtype=names, sep='\t').set_index(['chrom', 'start', 'end'])
    cscore2 = pd.read_csv(
        cscore2, usecols=[0,1,2,4], names=names.keys(),
        dtype=names, sep='\t').set_index(['chrom', 'start', 'end'])
    cscore = pd.merge(cscore1, cscore2, left_index=True, right_index=True).dropna()

    # Cscore sign is arbitrarily but we expect the same compartments to be
    # positively correlated without explicity knowing which is A or B.
    rho, p = pearsonr(cscore['cscore_x'], cscore['cscore_y'])
    if rho < 0:
        cscore['cscore_y'] *= -1
    cscore['switch'] = (cscore['cscore_x'] - cscore['cscore_y']).abs()
    cscore = cscore.reset_index()
    cscore['name'] = '.'
    cscore[['chrom', 'start', 'end', 'name', 'switch']].to_csv(
        sys.stdout, header=False, index=False, sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=computeSwitchScore)
    parser.add_argument(
        'cscore1', help='Output BED file of Cscore.')
    parser.add_argument(
        'cscore2', help='Output BED file of Cscore.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
