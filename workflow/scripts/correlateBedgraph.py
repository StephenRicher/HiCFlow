#!/usr/bin/env python3

""" Correlate 2 rescaled bedgraphs, in JSON format, with same window size. """

import sys
import json
import logging
import argparse
import pandas as pd
from scipy import stats
from typing import List
from utilities import setDefaults
from collections import defaultdict
from itertools import combinations

__version__ = '1.0.0'


def runCorrelation(bedGraphs: List):
    print('bedGraph1', 'bedGraph2', 'window', 'r', 'p', sep='\t')
    for bedGraph1, bedGraph2 in combinations(bedGraphs, 2):
        correlateBedgraph(bedGraph1, bedGraph2)


def correlateBedgraph(bedGraph1: str, bedGraph2: str):
    """ Correlate rescaled bedgraphs """

    merged, window = mergeRescaled(bedGraph1, bedGraph2)
    if not merged.empty:
        r, p = stats.pearsonr(merged.iloc[:,0], merged.iloc[:,1])
        print(bedGraph1, bedGraph2, window, r, p, sep='\t')


def mergeRescaled(bedGraph1: str, bedGraph2: str):
    """ Return bedgraph pair as merged pandas """
    bed1, window1 = processRescaled(bedGraph1)
    bed2, window2 = processRescaled(bedGraph2)
    if window1 != window2:
        logging.warning(
            f'Window sizes of {bedGraph1} and {bedGraph2} do not match.')
        return pd.DataFrame(), None
    else:
        return pd.concat([bed1, bed2], axis=1).dropna(), window1


def processRescaled(file: str):
    """ Process JSON format output of rescale bedgraph """

    with open(file) as fh:
        data = json.load(fh)
    df = (pd.DataFrame.from_dict(data['data'])
     .reset_index()
     .melt(id_vars='index', value_vars=data['data'],
           var_name='chromosome', value_name=file)
     .set_index(['chromosome', 'index']))
    df.index.names = ['chromosome', f'window-{data["window"]}']

    return df, data['window']


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'bedGraphs', nargs='*',
        help='Rescaled bedgraph files of equal window size, '
             'output out of rescaledBedgraph.py.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(runCorrelation(**vars(args)))
