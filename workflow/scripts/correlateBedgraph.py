#!/usr/bin/env python3

""" Correlate 2 rescaled bedgraphs, in JSON format, with same window size. """

import sys
import json
import logging
import argparse
import pandas as pd
import seaborn as sns
from scipy import stats
from typing import List
import statsmodels.api as sm
import matplotlib.pyplot as plt
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
        print(bedGraph1, bedGraph2, merged.attrs['window'], testStat, p, sep='\t')


def compareBedGraph(
        booleanBedGraph: str, bedGraph: str, test: str, normalityPlot: str):

    merged = mergeRescaled(booleanBedGraph, bedGraph)
    if merged.attrs['mode1'] != 'binary':
        raise FormatError(
            f'{booleanBedGraph} is not a boolean JSON bedGraph.')

    merged = merged.pivot(columns=booleanBedGraph)
    falseValues = merged[(bedGraph, False)].dropna()
    trueValues = merged[(bedGraph, True)].dropna()
    # Perform chi-square if other bedGraph is binary
    if merged.attrs['mode2'] == 'binary':
        obs = merged.apply(pd.Series.value_counts)
        testStat, p, dof, ex = stats.chi2_contingency(obs)
        if (ex.min().min()) < 5 or (obs.min().min() < 5):
            logging.error(
                'Atleast 1 of observed or expected values less than 5. '
                'Test result may be invalid.')
        test = 'chisquare'
    else:
        if normalityPlot is not None:
            plotNormalityCheck(trueValues, falseValues, normalityPlot)
        if test == 'mannwhitneyu':
            testStat, p = stats.mannwhitneyu(
                trueValues, falseValues)
        else:
            testStat, p = stats.ttest_ind(
                trueValues, falseValues, equal_var=True)

    # Write summary stats to stderr
    merged.columns = merged.columns.droplevel(0)
    summary = merged.describe().to_dict()

    testResult = {
        'groupData':     booleanBedGraph,
        'valueData':     bedGraph,
        'window':        merged.attrs['window'],
        'test':          test,
        'testStatistic': testStat,
        'p':             p,
        'summary':       summary}
    json.dump(testResult, sys.stdout, indent=4)


def plotNormalityCheck(trueValues, falseValues, out):
    """ Plot QQ-plot and histogram with Anderson darling test """

    # Anderson Darling Test for True values
    ad2, p = sm.stats.diagnostic.normal_ad(trueValues)
    # Probability plot of True group values
    ax1 = plt.subplot(221)
    sm.qqplot(trueValues, line='s', ax=ax1)
    ax1.set_title('QQ-plot')
    # Histgram of True group values
    ax2 = plt.subplot(222)
    sns.histplot(trueValues, kde=True)
    plt.text(1.1, 0.5, f'Group: True\n(ad2: {round(ad2, 1)}, p: {round(p, 3)})',
             transform=ax2.transAxes, rotation=270, fontsize='medium',
             verticalalignment='center', horizontalalignment='center')
    ax2.set_title('Histogram')
    ax2.set_xlabel('')

    # Anderson Darling Test for True values
    ad2, p = sm.stats.diagnostic.normal_ad(trueValues)
    # Probability plot of False group values
    ax3 = plt.subplot(223)
    sm.qqplot(falseValues, line='s', ax=ax3)
    ax3.set_title('')
    # Histgram of False group values
    ax4 = plt.subplot(224)
    sns.histplot(falseValues, kde=True)
    plt.text(1.1, 0.5, f'Group: True\n(ad2: {round(ad2, 1)}, p: {round(p, 3)})',
             transform=ax4.transAxes, rotation=270, fontsize='medium',
             verticalalignment='center', horizontalalignment='center')
    ax4.set_xlabel('')
    plt.tight_layout()
    plt.savefig(out, dpi=300)


def mergeRescaled(bedGraph1: str, bedGraph2: str):
    """ Return bedgraph pair as merged pandas """
    bed1 = processRescaled(bedGraph1)
    bed2 = processRescaled(bedGraph2)
    if bed1.attrs['window'] != bed2.attrs['window']:
        raise FormatError(
            f'Window sizes of {bedGraph1} and {bedGraph2} do not match.')
    else:
        merged = pd.concat([bed1, bed2], axis=1).dropna()
        merged.attrs['window'] = bed1.attrs['window']
        merged.attrs['mode1'] = bed1.attrs['mode']
        merged.attrs['mode2'] = bed2.attrs['mode']
        return merged


def processRescaled(file: str):
    """ Process JSON format output of rescale bedgraph """

    with open(file) as fh:
        data = json.load(fh)
    df = (pd.DataFrame.from_dict(data['data'])
         .reset_index()
         .melt(id_vars='index', value_vars=data['data'],
               var_name='chromosome', value_name=file)
         .set_index(['chromosome', 'index'])
         .dropna())
    df.index.names = ['chromosome', f'window-{data["window"]}']

    for attribute, value in data.items():
        if attribute == 'data':
            continue
        df.attrs[attribute] = value

    return df


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
        description=compareBedGraph.__doc__,
        help='Perform t-test or Mann Whitney-U on rescaled JSON bedgraph.',
        epilog=parser.epilog, parents=[mainParent])
    compare.add_argument(
        'booleanBedGraph', help='Rescaled boolean bedGraph to define groups.')
    compare.add_argument(
        'bedGraph', help='Rescaled bedGraph.')
    compare.add_argument(
        '--test', default='mannwhitneyu', choices=['ttest', 'mannwhitneyu'],
        help='Statistical test to apply (default: %(default)s)')
    compare.add_argument(
        '--normalityPlot',
        help='Write normality assessement plot to path (default: %(default)s)')
    compare.set_defaults(function=compareBedGraph)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
