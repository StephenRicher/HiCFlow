#!/usr/bin/env python3

""" Read loop domains and compute mean score at region """

import sys
import logging
import argparse
import fileinput
from scipy import stats
from typing import List
from utilities import setDefaults, createMainParent, readHomer


__version__ = '1.0.0'


def meanLoops(matrix: str, loops: List, binSize: int, absolute: bool, Z: bool):

    mat = readHomer(matrix, binSize, sparse=True)
    if absolute:
        mat['score'] = mat['score'].abs()
    with fileinput.input(loops) as fh:
        for line in fh:
            chr1, start1, end1, chr2, start2, end2 = parseLoops(line)
            if mat.attrs['chrom'] != chr1:
                continue
            overlap = getOverlapping(mat, start1, end1, start2, end2)
            if overlap.empty:
                continue
            score = overlap['score'].mean()
            if Z:
                # Round loopsize to nearest binSize
                loopSize = roundBin(start2 - start1, mat.attrs['binSize'])
                if loopSize == 0:
                    continue
                # Get all scores corresponding to same interaction distance
                regions = mat.loc[mat['start2'] - mat['start'] == loopSize, 'score']
                std = regions.std(ddof=0)
                mean = regions.mean()
                # Compute Z transform
                score = (score - mean) / std
            print(chr1, start1, end1, chr2, start2, end2, score, sep='\t')


def roundBin(x, base):
    return base * round(x/base)

def getOverlapping(mat, start1, end1, start2, end2):
    """ Return matrix interactions overlapping the interval pairs """

    # Loop1 overlaps if start or end is between matrix start/end interval
    endMat = mat['start'] + mat.attrs['binSize']
    start1LoopOverlap = (int(start1) >= mat['start']) & (int(start1) < endMat)
    end1LoopOverlap = (int(end1) >= mat['start']) & (int(end1) < endMat)
    loop1Overlap = start1LoopOverlap | end1LoopOverlap

    # Same for loop2
    endMat2 = mat['start2'] + mat.attrs['binSize']
    start2LoopOverlap = (int(start2) >= mat['start2']) & (int(start2) < endMat2)
    end2LoopOverlap = (int(end2) >= mat['start2']) & (int(end2) < endMat2)
    loop2Overlap = start2LoopOverlap | end2LoopOverlap

    return mat.loc[loop1Overlap & loop2Overlap]


def parseLoops(line):
    """ Reat bedgraph entry and convert types as appropriate """
    chr1, start1, end1, chr2, start2, end2 = line.split()[:6]
    return chr1, int(start1), int(end1), chr2, int(start2), int(end2)


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
             'computing mean. May be Appropriate for logFC '
             'comparison matrices. (default: %(default)s)')
    parser.add_argument(
        '--Z', action='store_true',
        help='Convert loop scores to Z scores. (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--binSize', required=True,
        type=int, help='Bin size for matrix.')
    parser.set_defaults(function=meanLoops)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
