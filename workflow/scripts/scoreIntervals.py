#!/usr/bin/env python3

""" Summarise bedgraph score per interval region """

import sys
import argparse
from collections import defaultdict
from bedgraphUtils import readBed
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def scoreIntervals(bedGraph: str, bed: str, buffer: int):
    bedgraph = readBed(bedGraph, filetype='bedgraph')
    records = readBed(bed, buffer)
    scoredRegions = defaultdict(float)
    for chrom, beds in records.items():
        for bed in beds:
            score = 0
            validRanges, remove = getValidRanges(bed, bedgraph[chrom])
            try:
                del bedgraph[chrom][:remove + 1]
            except TypeError:
                pass
            if not validRanges:
                continue
            bedInterval = set(bed.interval)
            for validRange in validRanges:
                # Detect base overlap between bedgraph interval and each region
                overlap = getOverlap(bedInterval, validRange.interval)
                score += validRange.normScore * overlap
            print(bed.chrom, bed.start, bed.end, bed.name, score, sep='\t')


def getOverlap(set1, range1):
    return len(set1.intersection(range1))


def getValidRanges(record, recordList):
    """ Return BED objects that overlap sorted list of BED objects.
        Also return upper index that no longer needs to be checked
        in next runs. """

    ranges = []
    minInterval = record.start
    maxInterval = record.end
    remove = None
    for i, bed in enumerate(recordList):
        if minInterval > bed.end:
            remove = i
        elif maxInterval < bed.start:
            break
        else:
            ranges.append(bed)

    return ranges, remove


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument(
        'bedGraph', help='BedGraph/BED interval file to rescale.')
    parser.add_argument(
        'bed', metavar='BED', help='BED file indicating regions to process.')
    parser.add_argument(
        '--buffer', type=int, default=0,
        help='Extend BED regions by +/- this value (default: %(default)s)')
    parser.set_defaults(function=scoreIntervals)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
