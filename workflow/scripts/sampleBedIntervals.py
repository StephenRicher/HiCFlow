#!/usr/bin/env python3

""" Randomly create BED intervals from referenceBED based on
    length in sampleBED  """

import sys
import random
import argparse
from typing import List
from collections import defaultdict
from utilities import setDefaults, createMainParent
from bedgraphUtils import splitPos, readRegions


__version__ = '1.0.0'


def sampleIntervals(referenceBed: str, sampleBed: str, nSamples: int):
    regions = readRegions(referenceBed)
    lengths = readLengths(sampleBed)
    lengths = randomiseLengths(lengths, nSamples)
    for length in lengths:
        chrom, start, end = getRandomPos(regions, length)
        print(chrom, start, end, sep='\t')


def randomiseLengths(lengths: List, n: int = None):
    """ Return random sample of 'n' lengths. Each length
        sampled equally such that if n = len(lengths) all
        values of length will be returned. """
    if n is None:
        n = len(lengths)
    lengthsCopy = []
    randomisedLengths = []
    for rep in range(n):
        if not randomisedLengths:
            lengthsCopy = lengths.copy()
            random.shuffle(lengthsCopy)
        randomisedLengths.append(lengthsCopy.pop())
    return randomisedLengths


def readLengths(bed):
    """ Return list of region length from BED input """
    lengths = []
    with open(bed) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chrom, start, end = splitPos(line)
            lengths.append(end - start)
    return lengths


def getRandomPos(regions, length=1, maxAttempts=100, _attempts=0):
    """ Extract random genomic start, end coordinates that fully
        overlap an interval. Repeat up to maxAttempts. """
    # Get total bases in intervals
    totalLength = 0
    for chrom, intervals in regions.items():
        for interval in intervals:
            totalLength += len(interval)
    # Generate random start index
    pos = random.randint(0, totalLength - 1)
    # Extract start position using index
    currentTotal = 0
    for chrom, intervals, in regions.items():
        for interval in intervals:
            if pos < currentTotal + len(interval):
                start = interval[(pos - currentTotal)]
                # If start + length not in same interval then repeat selection
                if pos + length >= currentTotal + len(interval):
                    if _attempts > maxAttempts:
                        return None
                    _attempts += 1
                    coords = getRandomPos(
                        regions, length, maxAttempts, _attempts)
                    if coords is None:
                        return None
                else:
                    coords = (chrom, start, start + length)
                return coords
            currentTotal += len(interval)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument(
        'referenceBed',
        help='BED file containing intervals within '
             'which regions will be sampled .')
    parser.add_argument(
        'sampleBed',
        help='BED file to extract interval lengths of sample.')
    parser.add_argument(
        '--nSamples', type=int,
        help='Number of intervals to sample (default: same as sampleBED).')
    parser.set_defaults(function=sampleIntervals)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
