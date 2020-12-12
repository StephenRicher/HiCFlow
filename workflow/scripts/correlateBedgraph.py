#!/usr/bin/env python3

""" Correlate 2 rescaled bedgraphs, in JSON format, with same window size. """

import sys
import json
import argparse
import pandas as pd
from scipy import stats
from typing import List
from utilities import setDefaults
from collections import defaultdict

__version__ = '1.0.0'


def correlateBedgraph(bedGraphs: List):
    """ Correlate rescaled bedgraphs """

    merged = pd.concat(
        [processRescaled(file) for file in bedGraphs], axis=1).dropna()
    assert merged.index.names != [None, None], f'Window sizes do not match'
    r, p = stats.pearsonr(merged.iloc[:,0], merged.iloc[:,1])
    sys.stdout.write(f'r\tp\n{r}\t{p}\n')


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

    return df


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'bedGraphs', nargs=2,
        help='Rescaled bedgraph files of equal window size, '
             'output out of rescaledBedgraph.py.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(correlateBedgraph(**vars(args)))
