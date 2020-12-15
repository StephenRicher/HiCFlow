#!/usr/bin/env python3

""" Correlate 2 rescaled bedgraphs, in JSON format, with same window size. """

import sys
import json
import logging
import argparse
import pandas as pd
from scipy import stats
from typing import List
from collections import defaultdict
from itertools import combinations
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'

def runCorrelation(bedGraphs: List):
    print('bedGraph1', 'bedGraph2', 'window', 'r', 'p', sep='\t')
    for bedGraph1, bedGraph2 in combinations(bedGraphs, 2):
        correlateBedgraph(bedGraph1, bedGraph2)


def correlateBedgraph(bedGraph1: str, bedGraph2: str):
    """ Correlate rescaled bedgraphs """

    merged = mergeRescaled(bedGraph1, bedGraph2)
    if not merged.empty:
        testStat, p = stats.pearsonr(merged.iloc[:,0], merged.iloc[:,1])
        print(bedGraph1, bedGraph2, merged['window'], testStat, p, sep='\t')


def runCompare(booleanBedGraph: str, bedGraphs: List, test: str):
    print('booleanBedGraph', 'bedGraph', 'window', 'test', 'testStatistic',
          'p', 'medianTrue', 'medianFalse', 'meanTrue', 'meanFalse',
          'countTrue', 'countFalse', 'sumTrue', 'sumFalse', sep='\t')
    for bedGraph in bedGraphs:
        compareBedGraph(booleanBedGraph, bedGraph, test)


def compareBedGraph(booleanBedGraph: str, bedGraph: str, test: str):

    merged, meta = mergeRescaled(booleanBedGraph, bedGraph)
    if not meta['binary1']:
        raise FormatError(
            f'{booleanBedGraph} is not a boolean JSON bedGraph.')

    merged = merged.pivot(columns=booleanBedGraph)
    falseValues = merged[(bedGraph, False)].dropna()
    trueValues = merged[(bedGraph, True)].dropna()
    # Perform chi-square if other bedGraph is binary
    if meta['binary2']:
        testStat, p, dof, ex = stats.chi2_contingency(
            merged.apply(pd.Series.value_counts))
        test = 'chisquare'
    else:
        if test == 'mannwhitneyu':
            testStat, p = stats.mannwhitneyu(
                trueValues, falseValues)
        else:
            testStat, p = stats.ttest_ind(
                trueValues, falseValues, equal_var=True)
    print(booleanBedGraph, bedGraph, meta['window'], test, testStat, p,
          trueValues.median(), falseValues.median(),
          trueValues.mean(), falseValues.mean(),
          trueValues.count(), falseValues.count(),
          trueValues.sum(), falseValues.sum(), sep='\t')


def mergeRescaled(bedGraph1: str, bedGraph2: str):
    """ Return bedgraph pair as merged pandas """
    bed1, meta1 = processRescaled(bedGraph1)
    bed2, meta2 = processRescaled(bedGraph2)
    if meta1['window'] != meta2['window']:
        raise FormatError(
            f'Window sizes of {bedGraph1} and {bedGraph2} do not match.')
    else:
        meta = {}
        merged = pd.concat([bed1, bed2], axis=1).dropna()
        meta['window'] = meta1['window']
        meta['binary1'] = meta1['binary']
        meta['binary2'] = meta2['binary']
        return merged, meta


def processRescaled(file: str):
    """ Process JSON format output of rescale bedgraph """
    meta = {}
    with open(file) as fh:
        data = json.load(fh)
    df = (pd.DataFrame.from_dict(data['data'])
     .reset_index()
     .melt(id_vars='index', value_vars=data['data'],
           var_name='chromosome', value_name=file)
     .set_index(['chromosome', 'index']))
    df.index.names = ['chromosome', f'window-{data["window"]}']
    meta['window'] = data['window']
    meta['binary'] = data['binary']

    return df, meta


class FormatError(Exception):
    pass


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    subparser = parser.add_subparsers(
        title='required commands', description='',
        metavar='Commands', help='Description:')
    correlate = subparser.add_parser(
        'correlate',
        description=correlateBedgraph.__doc__,
        help='Pairise correlation of bedgraphs.',
        epilog=parser.epilog, parents=[mainParent])
    correlate.add_argument(
        'bedGraphs', nargs='*',
        help='Rescaled bedgraph files of equal window size, '
             'output out of rescaledBedgraph.py.')
    correlate.set_defaults(function=runCorrelation)

    compare = subparser.add_parser(
        'compare',
        description=runCompare.__doc__,
        help='Perform t-test or Mann Whitney-U on rescaled JSON bedgraph.',
        epilog=parser.epilog, parents=[mainParent])
    compare.add_argument(
        'booleanBedGraph', help='Rescaled boolean bedGraph to define groups.')
    compare.add_argument(
        'bedGraphs', nargs='*', help='Rescaled bedGraph.')
    compare.add_argument(
        '--test', default='mannwhitneyu', choices=['ttest', 'mannwhitneyu'],
        help='Statistical test to apply.')
    compare.set_defaults(function=runCompare)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
