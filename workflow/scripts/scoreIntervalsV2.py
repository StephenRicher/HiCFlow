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
    bedgraph = sortBedgraph(bedgraph)
    print(bedgraph)
    regions = readBed(bed)
    scoredRegions = defaultdict(float)
    for chrom, intervals in regions.items():
        for name, interval in intervals:
            print(chrom, interval, list(bedgraph[chrom].keys()))
            validRanges = getValidRanges(interval, list(bedgraph[chrom].keys()))
            print(validRanges)
            for validRange in validRanges:
                # Detect base overlap between bedgraph interval and each region
                overlap = getOverlap(validRange, interval)
                key = f'{chrom} {min(interval)} {max(interval)+1} {name}'
                score = bedgraph[chrom][validRange] * len(overlap)
                scoredRegions[key] += score
    for region, score in scoredRegions.items():
        chrom, start, end, name = region.split()
        print(chrom, start, end, name, score, sep='\t')


def getOverlap(range1, range2):
    return range(max(min(range1), min(range2)), min(max(range1), max(range2))+1)


def getValidRanges(interval, rangeList):
    """ Return range objects that overlap the interval.
        Must provided sorted list of ranges. """

    ranges = []
    for i, rangeObj in enumerate(rangeList):
        if interval.start in rangeObj or interval.stop in rangeObj:
            ranges.append(rangeObj)
        if rangeObj.start > interval.stop:
            break
    return ranges


def readBedgraph(file):
    """ Construct bedgraph dictionary per chromosome with interval ranges """
    bedgraph = defaultdict(dict)
    with open(file) as fh:
        for line in fh:
            chrom, start, end, score = line.split()
            regionLength = int(end) - int(start)
            bedgraph[chrom][range(int(start),int(end))] = float(score) / regionLength
    return bedgraph


def sortBedgraph(bedgraph):
    """ Sort ranges in bedgraph dictionary per chromosome """
    sortedRegions = defaultdict(dict)
    for chrom, intervals in bedgraph.items():
        sortedIntervals = sorted(list(intervals), key=lambda r: r.start)
        for interval in sortedIntervals:
            sortedRegions[chrom][interval] = bedgraph[chrom][interval]
    return sortedRegions


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


def splitPos(line):
    chrom, start, end = line.split()[:3]
    return chrom, int(start), int(end)


def splitName(line):
    chrom, start, end, name = line.split()[:4]
    return chrom, int(start), int(end), name


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
