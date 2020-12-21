#!/usr/bin/env python3

""" Summarise bedgraph score per interval region """

import sys
import argparse
from collections import defaultdict
from bedgraphUtils import Bed, Bedgraph
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
            if remove is not None:
                del bedgraph[chrom][:remove+1]
            if not validRanges:
                continue
            for validRange in validRanges:
                # Detect base overlap between bedgraph interval and each region
                overlap = getOverlap(validRange.interval, bed.interval)
                score += validRange.normScore * len(overlap)
            print(bed.chrom, bed.start, bed.end,
                  bed.name, score, sep='\t', flush=True)


def getOverlap(range1, range2):
    return range(max(min(range1), min(range2)), min(max(range1), max(range2))+1)


def getValidRanges(record, recordList):
    """ Return BED objects that overlap sorted list of BED objects """

    ranges = []
    minInterval = record.start
    maxInterval = record.end
    remove = None
    for i, bed in enumerate(recordList):
        if minInterval > bed.end:
            remove = i
        elif minInterval in bed.interval or maxInterval in bed.interval:
            ranges.append(bed)
        elif bed.start > maxInterval:
            break
    return ranges, remove


def readBed(file, buffer=0, filetype='bed'):
    """ Construct sorted dictionary of Bed/Bedgraph objects """
    assert filetype in ['bed', 'bedgraph']
    records = defaultdict(list)
    with open(file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if filetype == 'bed':
                bed = Bed(line, buffer)
            else:
                bed = Bedgraph(line)
            records[bed.chrom].append(bed)
    # Sort per-chromosome ranges by start
    for chrom, bed in records.items():
        records[chrom] = sorted(bed, key=lambda r: r.start)
    return records


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
