#!/usr/bin/env python3

""" Summarise bedgraph score per interval region """

import sys
import argparse
from collections import defaultdict
from utilities import setDefaults, createMainParent
from bedgraphUtils import splitScore, splitPos, splitName


__version__ = '1.0.0'


def scoreIntervals(bedGraph: str, bed: str):
    bedgraph = readBedgraph(bedGraph)
    regions = readBed(bed)
    scoredRegions = defaultdict(float)
    for chrom, intervals in regions.items():
        for name, interval in intervals:
            key = f'{chrom} {min(interval)} {max(interval)+1} {name}'
            validRanges = getValidRanges(interval, list(bedgraph[chrom].keys()))
            if not validRanges:
                continue
            for validRange in validRanges:
                # Detect base overlap between bedgraph interval and each region
                overlap = getOverlap(validRange, interval)
                score = bedgraph[chrom][validRange] * len(overlap)
                scoredRegions[key] += score
            print(chrom, min(interval), max(interval)+1,
                  name, scoredRegions[key], sep='\t')


def getOverlap(range1, range2):
    return range(max(min(range1), min(range2)), min(max(range1), max(range2))+1)


def getValidRanges(interval, rangeList):
    """ Return range objects that overlap the interval.
        Must provided sorted list of ranges. """

    ranges = []
    minInterval = min(interval)
    maxInterval = max(interval)
    for i, rangeObj in enumerate(rangeList):
        if minInterval in rangeObj or maxInterval in rangeObj:
            ranges.append(rangeObj)
        if min(rangeObj) > maxInterval:
            break
    return ranges


def readBedgraph(file):
    """ Construct bedgraph dictionary per chromosome with interval ranges """
    bedgraph = defaultdict(dict)
    with open(file) as fh:
        for line in fh:
            chrom, start, end, score = splitScore(line, filetype='bedgraph')
            regionLength = end - start
            bedgraph[chrom][range(start, end)] = score / regionLength
    # Sort per-chromosome ranges by start
    for chrom, intervals in bedgraph.items():
        bedgraph[chrom] = dict(sorted(intervals.items(), key=lambda r: r[0].start))
    return bedgraph


def readBed(bed):
    """ Read bed file into dictionary structure """
    regions = defaultdict(list)
    with open(bed) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                chrom, start, end, name = splitName(line)
            except ValueError:
                chrom, start, end = splitPos(line)
                name = '.'
            regions[chrom].append((name, range(int(start), int(end))))
    return regions


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument(
        'bedGraph', help='BedGraph/BED interval file to rescale.')
    parser.add_argument(
        'bed', metavar='BED', help='BED file indicating regions to process.')
    parser.set_defaults(function=scoreIntervals)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
