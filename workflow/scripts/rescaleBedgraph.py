#!/usr/bin/env python3

""" Convert bedgraph intervals to new window size and write to JSON. """

import sys
import json
import logging
import argparse
from collections import defaultdict
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def rescaleBedgraph(
        bedGraph: str, chromSizes: str, window: int, regions: str, bed: bool,
        distanceTransform: bool, binary: bool, threshold: int, includeZero: bool):
    """ Convert intervals to new window size and write to JSON. """

    chromSizes = readChromSizes(chromSizes)
    scores = readBedgraph(bedGraph, window, bed, binary, threshold)
    regions = readRegions(regions)

    # Non-zero interval scores stored in dictionary
    rescaledBedgraph = {'window' :   window,
                        'binary' :   binary,
                        'threshold': threshold,
                        'data'   :   defaultdict(dict)}

    # Loop through windows (in order) per chromosome and perform distance
    # transformation as required.
    for chrom, size in chromSizes.items():
        # Reset prevScore for each chromosome
        prevScore = None
        for start in range(0, size, window):
            # Skip windows where start not in regions
            if regions and not validRegion(regions, chrom, start):
                prevScore = None
                continue
            try:
                score = scores[chrom][start]
            except KeyError:
                prevScore = 0
                if binary or includeZero:
                    score = 0
                else:
                    continue
            if binary:
                rescaledBedgraph['data'][chrom][start] = score > 0.5
            elif distanceTransform:
                try:
                    rescaledBedgraph['data'][chrom][start] = score - prevScore
                except TypeError:
                    pass # Skip windows where prevScore set to None
            else:
                rescaledBedgraph['data'][chrom][start] = score
            prevScore = score
    json.dump(rescaledBedgraph, sys.stdout)


def readRegions(bed):
    """ Read BED file and return dict of chromosomes and allowed intervals """
    regions = defaultdict(list)
    if bed is not None:
        with open(bed) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                chrom, start, end, score = splitBed(line)
                regions[chrom].append(range(start, end))
    return regions


def validRegion(regions, chrom, start):
    """ Return True if interval present in regions dict """
    for interval in regions[chrom]:
        if start in interval:
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
    except ValueError:
        chrom, start, end = line.split()[:3]
        score = 1
    return chrom, int(start), int(end), float(score)


def getWindow(pos, window):
    """ Return 0-based window start for a given position"""
    return (pos // window) * window


def readBedgraph(file, window, bed=False, binary=False, threshold=None):
    """ Read per-base scores into new window  """
    scores = defaultdict(dict)
    with open(file) as fh:
        for i, line in enumerate(fh):
            if bed: # Input is BED format
                # Check 1st line if score column included
                if i == 0 and noBedScore(line):
                    if not binary:
                        logging.error(
                            f'{file} has no score column, '
                             'Must activate "--binary" mode.')
                        sys.exit(1)
                    if threshold is not None:
                        logging.warning(
                            f'{file} has no score column, '
                             '--threshold may not be meaningul.')
                chrom, start, end, score = splitBed(line)
            else: # Input is bedgraph format
                chrom, start, end, score = splitBedgraph(line)
            # Set score to 1 if binary is set
            if binary:
                if threshold is not None:
                    score = score > threshold
                else:
                    score = True
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


def noBedScore(line):
    """ Return True if BED file contains no score column. """
    return len(line.split()) < 5


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=rescaleBedgraph)
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
        help='Perform differencing to remove series dependence. '
             'Ignored if --binary is set (default: %(default)s)')
    parser.add_argument(
        '--regions', help='BED file indicating regions to process.')
    parser.add_argument(
        '--bed', action='store_true',
        help='Treat input as BED format (default: %(default)s)')
    parser.add_argument(
        '--binary', action='store_true',
        help='Set boolean intervals e.g. presence/absence of a gene. If used, '
             '--includeZero is switched on (default: %(default)s)')
    parser.add_argument(
        '--threshold', type=float,
        help='Score threshold for determining binary intervals. '
             'Only appicable if --binary is set (default: %(default)s)')
    parser.add_argument(
        '--includeZero', action='store_true',
        help='Include windows with a 0 score (default: %(default)s)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
