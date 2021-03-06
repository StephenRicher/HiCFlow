#!/usr/bin/env python3

""" Read loop domains and compute mean score at region """

import sys
import logging
import argparse
import fileinput
import pandas as pd
from scipy import stats
from typing import List
from utilities import setDefaults, createMainParent, readHomer
try:
    from pandarallel import pandarallel
except ModuleNotFoundError:
    pass

__version__ = '1.0.0'


def sumLoops(matrix: str, loops: List, absolute: bool, distanceNorm: bool, threads: int):

    # Read as sparse since we are computing sum
    mat = readHomer(matrix, diagonal=False, upperOnly=True, sparse=True,
                    distanceNorm=distanceNorm, absolute=absolute)
    mat['end'] = mat['start'] + mat.attrs['binSize']
    mat['end2'] = mat['start2'] + mat.attrs['binSize']
    loops = pd.concat(
        readLoops(loop, chrom=mat.attrs['chrom'], cis=True) for loop in loops)
    # Remove interactions larger than maximum loop size
    maxLoopSize = abs(loops['end2'] -  loops['start1']).max()
    mat = mat.loc[abs(mat['start2'] - mat['start']) <= maxLoopSize]
    if (threads > 1) and ('pandarallel' in sys.modules):
        pandarallel.initialize(nb_workers=threads, verbose=0)
        loops['score'] = loops.parallel_apply(
            getOverlapping, args=([mat]), axis=1)
    else:
        loops['score'] = loops.apply(getOverlapping, args=([mat]), axis=1)
    loops.to_csv(sys.stdout, sep='\t', index=False, header=False)


def getOverlapping(loop, mat):
    """ Return matrix interactions overlapping the interval pairs """

    # Loop1 overlaps if start or end is between matrix start/end interval
    start1LoopOverlap = (loop['start1'] >= mat['start']) & (loop['start1'] < mat['end'])
    end1LoopOverlap = (loop['end1'] >= mat['start']) & (loop['end1'] < mat['end'])
    loop1Overlap = start1LoopOverlap | end1LoopOverlap
    # Same for loop2
    start2LoopOverlap = (loop['start2'] >= mat['start2']) & (loop['start2'] < mat['end2'])
    end2LoopOverlap = (loop['end2'] >= mat['start2']) & (loop['end2'] < mat['end2'])
    loop2Overlap = start2LoopOverlap | end2LoopOverlap
    return mat.loc[loop1Overlap & loop2Overlap, 'score'].sum()


def readLoops(file: str, chrom=None, cis=True):
    """ Read BED file into Pandas """
    columns = {0: 'chrom1', 1: 'start1', 2: 'end1',
               3: 'chrom2', 4: 'start2', 5: 'end2'}
    loops = pd.read_csv(file, sep='\t', comment='#', header=None)
    if len(loops.columns) > 6:
        loops = loops[[0,1,2,3,4,5]]
    loops = loops.rename(columns=columns)
    if cis:
        loops = loops.loc[loops['chrom1'] == loops['chrom2']]
    if chrom is not None:
        loops = loops.loc[loops['chrom1'] == chrom]
    # Reorder any entries where start2 < start1 (so they fall in upper matrix)
    revLoops = loops.loc[loops['start2'] < loops['start1']]
    revLoops.columns = ['chrom2', 'start2', 'end2', 'chrom1', 'start1', 'end1']
    loops = loops.loc[loops['start2'] >= loops['start1']].append(revLoops)
    return loops


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument(
        'matrix', help='Matrix file in HOMER format.')
    parser.add_argument(
        'loops', nargs='*', default=[],
        help='Loop domains in LINKS format (default: stdin)')
    parser.add_argument(
        '--absolute', action='store_true',
        help='Convert matrix score to absolute values before '
             'computing mean. Absolute values calculated prior to '
             'distance normalisation. May be appropriate for logFC '
             'comparison matrices. (default: %(default)s)')
    parser.add_argument(
        '--distanceNorm', action='store_true',
        help='Normalise by interaction distance (obs/exp). Not '
             'recommended for comparison matrices (default: %(default)s)')
    parser.add_argument(
        '--threads', type=int, default=1,
        help='Threads for parallel processing (default: %(default)s)')
    parser.set_defaults(function=sumLoops)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
