#!/usr/bin/env python3

""" Randomly create BED intervals from referenceBED based on
    length in sampleBED  """

import sys
import random
import argparse
from itertools import repeat
from utilities import setDefaults, createMainParent
from bedgraphUtils import splitPos, readRegions


__version__ = '1.0.0'


def sampleIntervals(referenceBed: str, sampleBed: str, nRepeats: int):
    regions = readRegions(referenceBed)
    intervalLengths = readLengths(sampleBed)
    for i, lengths in enumerate(repeat(intervalLengths, nRepeats)):
        for length in lengths:
            chrom, start, end = getRandomPos(regions, length)
            print(chrom, start, end, i, sep='\t')


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
    """ Extract random genomic start/end coordinates that fully
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
        '--nRepeats', default=1, type=int,
        help='Number of repeat samples (default: %(default)s).')
    parser.set_defaults(function=sampleIntervals)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
