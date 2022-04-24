#!/usr/bin/env python3

""" Compute subtraction score at each looping domain """

import os
import sys
import argparse
import contextlib
import numpy as np
import pandas as pd
from typing import List
from pathlib import Path
from hicmatrix import HiCMatrix as hm
from sklearn.preprocessing import KBinsDiscretizer
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def scoreLoopDiff(loops: List, matrix: str, maxLineWidth: int,
                  nBins: int, interactOut: str, linksUp: str, linksDown: str):

    loops = readLoops(loops)
    matrix = hm.hiCMatrix(matrix)
    chrom = matrix.getChrNames()[0]
    # Remove non-specific loops
    loops = loops.loc[loops['chrom1'] == chrom]

    if loops.empty:
        for file in [interactOut, linksUp, linksDown]:
            Path(file).touch()
        return 0
        
    loops[['rawScore', 'direction']] = loops.apply(
        scoreLoops, axis=1, args=(matrix,), result_type='expand')
    loops = loops.dropna()

    nBins = min(nBins, len(loops))
    est = KBinsDiscretizer(n_bins=nBins, encode='ordinal', strategy='kmeans')
    if nBins > 1:
        loops['score'] = est.fit_transform(
            loops['rawScore'].to_numpy().reshape(-1,1)).astype(int)
    else:
        loops['score'] = 1
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

    writeLinks(loops.loc[(loops['direction'] == 1)], linksUp)
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


def scoreLoops(loop, matrix):
    indices1 = matrix.getRegionBinRange(loop['chrom1'], loop['start1'], loop['end1'])
    indices2 = matrix.getRegionBinRange(loop['chrom2'], loop['start2'], loop['end2'])
    try:
        score = matrix.matrix[indices1[0]-1:indices1[1]+1, indices2[0]-1:indices2[1]+1].toarray()
    except TypeError:
        # Loop out of range of matrix
        return np.nan, np.nan
    direction = -1 if score.mean() < 0 else 1
    score = np.absolute(score).mean()
    return score, direction


def readLoops(paths):
    names = ({
        'chrom1': str, 'start1': int, 'end1': int,
        'chrom2': str, 'start2': int, 'end2': int
    })
    loops = []
    for file in paths:
        loops.append(
            pd.read_csv(file, sep='\t', usecols=[0,1,2,3,4,5],
            names=names.keys(), dtype=names))
    loops = pd.concat(loops).drop_duplicates()
    loops = loops.loc[loops['chrom1'] == loops['chrom2']]
    return loops


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=scoreLoopDiff)
    parser.add_argument('loops', nargs='+', help='HiC loops file.')
    parser.add_argument(
        '--maxLineWidth', type=int, default=3,
        help='Maximum loop line width to plot for top '
             'scoring loop (default: %(default)s)')
    parser.add_argument(
        '--nBins', type=int, default=20,
        help='Number of clusters to group loops scores (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--matrix', help='HiC subtraction matrix in h5 format.')
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
