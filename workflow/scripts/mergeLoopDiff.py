#!/usr/bin/env python3

""" Merge output of scoreLoopDiff and scale scoring for plotting """

import os
import sys
import argparse
import numpy as np
import pandas as pd
from typing import List
from sklearn.preprocessing import KBinsDiscretizer
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def computeChangeScore(loops: List, maxLineWidth: int, nBins: int,
                       interactOut: str, linksUp: str, linksDown: str):
    loops = pd.concat([pd.read_pickle(loop) for loop in loops])
    est = KBinsDiscretizer(n_bins=nBins, encode='ordinal', strategy='kmeans')
    loops['score'] = est.fit_transform(
        loops['rawScore'].to_numpy().reshape(-1,1)).astype(int)
    loops['empty'] = '.'
    loops['color'] = 0
    loops['name'] = loops.index
    loops['size'] = (loops['start2'] - loops['start1']).abs()
    # Ensure start1 is less than start2 and otherwise swap
    s = loops['start2'] < loops['start1']
    loops.loc[s, ['start1','start2']] = loops.loc[s, ['start2','start1']].values
    loops.loc[s, ['end1','end2']] = loops.loc[s, ['end2','end1']].values
    # Ensure all integer
    loops[['start1', 'end1', 'start2', 'end2']] = (
        loops[['start1', 'end1', 'start2', 'end2']].astype(int))

    # Write full file to interact format
    writeInteract(loops, interactOut)

    # Rescale score to set max line width
    maxScore = (maxLineWidth * 2) ** 2
    loops['score'] = (loops['score'] / loops['score'].max()) * maxScore

    writeLinks(loops.loc[(loops['direction'] == 1)], linkUp)
    writeLinks(loops.loc[(loops['direction'] == -1)], linksDown)



def writeInteract(loops, out):
    """ Write loops file in Interact format for UCSC genome browser """
    names = ([
        'chrom1', 'start1', 'end2', 'name', 'score',
        'rawScore', 'empty', 'color',
        'chrom1', 'start1', 'end1', 'empty', 'empty',
        'chrom2', 'start2', 'end2', 'empty', 'empty'
    ])
    with contextlib.suppress(FileNotFoundError):
        os.remove(out)
    # Rescale score to 0 - 1000
    loops['score'] = loops['score'] / loops['score'].max() * 1000
    with open(out , 'a') as fh:
        fh.write(f'track type=interact useScore=on maxHeightPixels=200:100:50 '
                  'visibility=full\n')
        loops[names].to_csv(fh, sep='\t', header=False, index=False)


def writeLinks(loops, out):
    """ Write loops file in Links format for pyGenomeTracks """
    names = (['chrom1', 'start1', 'end2', 'chrom2', 'start2', 'end2', 'score'])
    loops[names].to_csv(out, sep='\t', header=False, index=False)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=computeChangeScore)
    parser.add_argument(
        'loops', nargs='+',
        help='Scored loop output files from scoreLoopDiff.')
    parser.add_argument(
        '--maxLineWidth', type=int, default=3,
        help='Maximum loop line width to plot for top '
             'scoring loop (default: %(default)s)')
    parser.add_argument(
        '--nBins', type=int, default=20,
        help='Number of clusters to group loops scores (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--interactOut', required=True,
        help='Path to Interact format output.')
    requiredNamed.add_argument(
        '--linksUp', required=True,
        help='Path to Links format positive change loop data.')
    requiredNamed.add_argument(
        '--linksDown', required=True,
        help='Path to Links format negative change loop data.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
