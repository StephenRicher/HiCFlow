#!/usr/bin/env python3

""" Plot Z-score distributed scored BED regions. Use 1 BED as a reference """


import sys
import argparse
import pandas as pd
import seaborn as sns
from typing import List
from scipy.stats import zscore
import matplotlib.pyplot as plt
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def plotBed(refBed: str, beds: List, out: str, dpi: int):
    infiles = []
    for rep, bed in enumerate(beds):
        infiles.append(processDat(bed, rep))
    infiles.append(processDat(refBed, 'reference'))
    meanScore = pd.concat(infiles)
    meanScore['z'] = zscore(meanScore['normScore'])

    fig, ax = plt.subplots()
    sns.histplot(meanScore['z'], kde=True, ax=ax)
    ax.axvline(meanScore.loc['reference','z'], color='red')
    ax.set_title(f'Z score normalised allelic differences.\nCpG islands '
                 f'(red line) vs random (n = {len(beds)})')
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)


def processDat(file: str, name: str):
    """ Read scored BED and summarise length normalised score """

    dtypes = {'chrom': 'str', 'start': 'int',
              'end': 'int', 'name': 'str',
              'score': 'float'}
    df = pd.read_csv(file, sep='\t', names=dtypes.keys(), dtype=dtypes)
    df['name'] = name
    df['normScore'] = (df['score']) / (df['end'] - df['start'])
    return df[['normScore','name']].groupby('name').mean()


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument(
        'refBed',
        help='Scored reference BED to compare against all others.')
    parser.add_argument(
        'beds', metavar='BED', nargs='+',
        help='Scored BED files to produce distribution.')
    parser.add_argument(
        '--dpi', type=int, default=300,
        help='DPI for figure (default: %(default)s)')
    required = parser.add_argument_group('required named arguments')
    required.add_argument(
        '--out', required=True, help='Output filename for plot.')
    parser.set_defaults(function=plotBed)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
