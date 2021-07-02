#!/usr/bin/env python3

""" Perform permutation testing by merging BAMs from both groups and randomly
    split (keeping the relative sample sizes approximately equal). Output is
    written to SUTM format and is repeated a determined number of times. """

import os
import sys
import argparse
import numpy as np
import pandas as pd
from typing import List
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def permutationTest(trueLogFC: str, allShadow: List, rawOut: str, binSize: int, chrom: str):

    trueLogFC = pd.read_pickle(trueLogFC)
    trueLogFC['shadow'] = 0

    # To save memory each file could be opened and read in chunksizes of (N x nShadow)
    # Process N start positions at a time
    # See https://stackoverflow.com/questions/59983073/how-to-load-pickle-file-in-chunks
    # The pickle would need to be saved in chunks within shadowSUTM and then iteratively read
    # Could write a seperate wscript to read pickles and merge 

    # Read and concanate all shadows, correctly adjust rep numbers
    nShadow = 0
    allData = []
    for shadow in allShadow:
        shadow = pd.read_pickle(shadow)
        shadow['shadow'] = shadow['shadow'] + nShadow
        nShadow = shadow['shadow'].max()
        allData.append(shadow)
    allData = pd.concat(allData + [trueLogFC])

    permutation = (allData.groupby('start1').apply(runPermute)
        .reset_index().rename({0: 'p'}, axis=1))

    permutation['name'] = '.'
    permutation['chrom'] = chrom
    permutation['end'] = permutation['start1'] + binSize
    columns = ['chrom', 'start1', 'end', 'name', 'p']
    permutation[columns].to_csv(sys.stdout, index=False, header=False, sep='\t')

    if rawOut is not None:
        allData['chrom'] = chrom
        allData.to_csv(rawOut, header=True, sep='\t')


def runPermute(x):
    """ Compare shadow abs(score) against normal matrix """
    shadowScores = x.loc[x['shadow'] > 0, 'logFC']
    normalScore = float(x.loc[x['shadow'] == 0, 'logFC'])
    totalAbove = (shadowScores >= normalScore).sum()
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
