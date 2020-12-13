#!/usr/bin/env python3

""" Convert bedgraph intervals to new window size and write to JSON. """

import sys
import json
import argparse
from utilities import setDefaults
from collections import defaultdict

__version__ = '1.0.0'


def rescaleBedgraph(
        bedGraph: str, chromSizes: str, window: int,
        distanceTransform: bool):
    """ Convert intervals to new window size and write to JSON. """

    chromSizes = readChromSizes(chromSizes)
    scores = readBedgraph(bedGraph, window)
    # Non-zero itnervals stored in dictionary
    rescaledBedgraph = {'window' : window,
                        'data'   : defaultdict(dict)}
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



def splitBedgraph(line):
    """ Split bedgraph columns and set type """
    chrom, start, end, score = line.split()
    return chrom, int(start), int(end), float(score)


def getWindow(pos, window):
    """ Return 0-based window start for a given position"""
    return (pos // window) * window


def readBedgraph(file, window):
    """ Read per-base scores into new window  """
    scores = defaultdict(dict)
    with open(file) as fh:
        for line in fh:
            chrom, start, end, score = splitBedgraph(line)
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
            chrom, size = line.split()
            chromSizes[chrom] = int(size)
    return chromSizes


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'bedGraph', help='BedGraph interval file to rescale.')
    parser.add_argument(
        'chromSizes', help='Chromosome sizes file.')
    parser.add_argument(
        '--window', type=int, default=100,
        help='Bedgraph interval window to rescale to (default: %(default)s)')
    parser.add_argument(
        '--distanceTransform', action='store_true',
        help='Perform differencing to remove series dependence.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(rescaleBedgraph(**vars(args)))
