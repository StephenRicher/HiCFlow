#!/usr/bin/env python3

""" Randomly create BED intervals from referenceBED based on
    length in sampleBED  """

import sys
import random
import logging
import argparse
from bedgraphUtils import readBed
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def sampleIntervals(bed: str, nRepeats: int, length: int, seed: float):
    random.seed(seed)
    regions = readBed(bed, filetype='bed')
    totalLength = getTotalLength(regions)
    for repeat in range(nRepeats):
        chrom, start, end = getRandomPos(regions, totalLength, length)
        print(chrom, start, end, repeat, sep='\t')


def getRandomPos(regions, totalLength, length=1):
    """ Extract random genomic start/end coordinates
        from set of BED intervals. """

    # Generate random start index
    pos = random.randint(0, totalLength - 1)
    # Extract start position using index
    currentTotal = 0
    for chrom, intervals, in regions.items():
        for interval in intervals:
            if pos < currentTotal + interval.regionLength:
                mid = interval.interval[(pos - currentTotal)]
                a, b = splitNum(length)
                return chrom, mid - a, mid + b
            currentTotal += interval.regionLength
    return None


def getTotalLength(regions):
    totalLength = 0
    for chrom, intervals in regions.items():
        for interval in intervals:
            totalLength += interval.regionLength
    return totalLength


def splitNum(num):
    """ Split integer into 2 integers """
    if num % 2 == 0:
        splitNum = (int(num / 2), int(num / 2))
    else:
        splitNum = int(num // 2), int((num // 2) + 1)
    return splitNum


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument(
        'bed', metavar='BED',
        help='BED intervals within which regions will be sampled .')
    parser.add_argument(
        '--length', default=1, type=int,
        help='The length of the intervals to generate.')
    parser.add_argument(
        '--nRepeats', default=100_000, type=int,
        help='Number of intervals to generate samples (default: %(default)s).')
    parser.add_argument(
        '--seed', default=None, type=float,
        help='Seed for random number generation (default: %(default)s)')
    parser.set_defaults(function=sampleIntervals)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
