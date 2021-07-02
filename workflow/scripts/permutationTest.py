#!/usr/bin/env python3

""" Perform permutation testing by merging BAMs from both groups and randomly
    split (keeping the relative sample sizes approximately equal). Output is
    written to SUTM format and is repeated a determined number of times. """

import os
import sys
import pickle
import argparse
import numpy as np
import pandas as pd
from typing import List
from contextlib import ExitStack
from statsmodels.stats.multitest import fdrcorrection
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def permutationTest(trueLogFC: str, allShadow: List, rawOut: str, binSize: int, chrom: str):

    allFiles = [trueLogFC] + allShadow
    allPermute = []
    with ExitStack() as stack:
        files = [stack.enter_context(open(file, 'rb')) for file in allFiles]
        while True:
            try:
                allData = pd.concat([pickle.load(fh) for fh in files])
                permutation = (allData.groupby('start1').apply(runPermute)
                        .reset_index().rename({0: 'p'}, axis=1))
                allPermute.append(permutation)
            except EOFError:
                break
    allPermute = pd.concat(allPermute)

    allPermute ['fdr'] =  fdrcorrection(allPermute['p'])[1]
    allPermute['chrom'] = chrom
    allPermute['end'] = allPermute['start1'] + binSize
    columns = ['chrom', 'start1', 'end', 'p', 'fdr']
    allPermute[columns].to_csv(sys.stdout, index=False, header=False, sep='\t')


def runPermute(x):
    """ Compare shadow abs(score) against normal matrix """
    shadowScores = x.loc[x['shadow'] > 0, 'logFC']
    normalScore = float(x.loc[x['shadow'] == 0, 'logFC'])
    totalAbove = (shadowScores <= normalScore).sum()
    total = len(shadowScores)
    p = (totalAbove / total)
    if p == 0:
        p += (1 / total)
    return p


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'trueLogFC', help='Real adjusted logFC frequences (.pkl)')
    parser.add_argument(
        'allShadow', nargs='*', help='Shadow logFC frequences (.pkl)')
    parser.add_argument(
        '--chrom', help='Reference of correponding matrices.')
    parser.add_argument(
        '--binSize', type=int, help='Bin size of corresponding matrices')
    parser.add_argument(
        '--rawOut', help='Path for raw, unsummarised data.')
    parser.set_defaults(function=permutationTest)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
