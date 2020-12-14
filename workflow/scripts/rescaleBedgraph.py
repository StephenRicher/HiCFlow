#!/usr/bin/env python3

""" Convert bedgraph intervals to new window size and write to JSON. """

import sys
import json
import argparse
from utilities import setDefaults
from collections import defaultdict

__version__ = '1.0.0'


def rescaleBedgraph(
        bedGraph: str, chromSizes: str, window: int, regions: str, bed: bool,
        distanceTransform: bool, binary: bool, includeZero: bool):
    """ Convert intervals to new window size and write to JSON. """

    chromSizes = readChromSizes(chromSizes)
    scores = readBedgraph(bedGraph, window, bed, binary)
    regions = readRegions(regions)

    # Non-zero interval scores stored in dictionary
    rescaledBedgraph = {'window' : window,
                        'data'   : defaultdict(dict)}

    # Loop through windows (in order) per chromosome and perform distance
    # transformation as required.
    for chrom, size in chromSizes.items():
        # Reset prevScore for each chromosome
        try:
            del prevScore
        except NameError:
            pass
        for start in range(0, size, window):
            try:
                score = scores[chrom][start]
            except KeyError:
                prevScore = 0
                if includeZero:
                    score = 0
                else:
                    continue
            if not validRegion(regions, chrom, start, start + window):
                continue
            if distanceTransform:
                try:
                    rescaledBedgraph['data'][chrom][start] = score - prevScore
                except NameError:
                    # Skip first window where prevScore undefined
                    pass
            else:
                rescaledBedgraph['data'][chrom][start] = score
            prevScore = score
    json.dump(rescaledBedgraph, sys.stdout)


def readRegions(bed):
    """ Read BED file and return dict of chromosomes and allowed intervals """
    if bed is None:
        return None
    regions = defaultdict(list)
    with open(bed) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chrom, start, end, score = splitBed(line)
            regions[chrom].append(range(start,end))
    return regions


def validRegion(regions, chrom, start, end):
    """ Return True if interval present in regions dict """
    if regions is None:
        return True
    for interval in regions[chrom]:
        if (start in interval) and (end in interval):
            return True
    return False


def splitBedgraph(line):
    """ Split bedgraph columns and set type """
    chrom, start, end, score = line.split()
    return chrom, int(start), int(end), float(score)


def splitBed(line):
    """ Split BED columns and set type """
    try:
        chrom, start, end, name, score = line.split()[:5]
    except ValueError: # Set score to 0 for BED3 format
        chrom, start, end = line.split()[:3]
        score = 0
    return chrom, int(start), int(end), float(score)


def getWindow(pos, window):
    """ Return 0-based window start for a given position"""
    return (pos // window) * window


def readBedgraph(file, window, bed=False, binary=False):
    """ Read per-base scores into new window  """
    scores = defaultdict(dict)
    with open(file) as fh:
        for line in fh:
            if bed: # Input is BED format
                chrom, start, end, score = splitBed(line)
            else: # Input is bedgraph format
                chrom, start, end, score = splitBedgraph(line)
            # Set score to 1 if binary is set
            if binary:
                score = 1
            regionLength = end - start
            for base in range(start, end):
                pos = getWindow(base, window)
                try:
                    scores[chrom][pos] += score / regionLength
                except KeyError:
                    scores[chrom][pos] = score / regionLength
    return scores


def readChromSizes(file):
    """ Read chromosome sizes to dict """
    chromSizes = {}
    with open(file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chrom, size = line.split()
            chromSizes[chrom] = int(size)
    return chromSizes


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'bedGraph',
        help='BedGraph interval file to rescale. '
             'For BED input must also set --bed.')
    parser.add_argument(
        'chromSizes', help='Chromosome sizes file.')
    parser.add_argument(
        '--window', type=int, default=100,
        help='Bedgraph interval window to rescale to (default: %(default)s)')
    parser.add_argument(
        '--distanceTransform', action='store_true',
        help='Perform differencing to remove series dependence.')
    parser.add_argument(
        '--regions', help='BED file indicating regions to process.')
    parser.add_argument(
        '--bed', action='store_true',
        help='Treat input as BED format (default: %(default)s)')
    parser.add_argument(
        '--binary', action='store_true',
        help='Set all scores 1. Useful for files that indicate boolean '
             'intervals e.g. presence/absence of a gene (default: %(default)s)')
    parser.add_argument(
        '--includeZero', action='store_true',
        help='Include windows with a 0 score (default: %(default)s)')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(rescaleBedgraph(**vars(args)))
