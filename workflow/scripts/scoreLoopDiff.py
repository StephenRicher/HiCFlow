#!/usr/bin/env python3

""" Compute subtraction score at each looping domain """

import os
import sys
import argparse
import numpy as np
import pandas as pd
from typing import List
from hicmatrix import HiCMatrix as hm
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def scoreLoopDiff(loops: List, matrix: str, out: str):

    loops = readLoops(loops)
    matrix = hm.hiCMatrix(matrix)
    chrom = hic.getChrNames()[0]
    # Remove non-specific loops
    loops = loops.loc[loops['chrom1'] == chrom]
    loops[['score', 'direction']] = loops.apply(
        scoreLoops, axis=1, args=(matrix,), result_type='expand')
    loops.dropna().to_pickle(out)


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
    loops = pd.concat(loops).unique()
    loops = loops.loc[loops['chrom1'] == loops['chrom2']]
    return loops


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=scoreLoopDiff)
    parser.add_argument('loops', nargs='+', help='HiC loops file.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--matrix', help='HiC subtraction matrix in h5 format.')
    requiredNamed.add_argument(
        '--out', help='Path to write pickled loop data.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
