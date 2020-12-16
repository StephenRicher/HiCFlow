#!/usr/bin/env python3

""" Convert bedgraph intervals to new window size and write to JSON. """

import sys
import json
import logging
import argparse
from collections import defaultdict
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def rescaleBinCount(bedGraph: str, chromSizes: str, window: int,
                    regions: str, threshold: float = None):

    scores = readBedgraph(bedGraph, window, threshold)
    mode = 'binary' if threshold is not None else 'count'
    template = makeTemplate(
        window, mode, regions, chromSizes,
        distanceTransform=False, includeZero=True, threshold=threshold)

    json.dump(makeRescaled(scores, template), sys.stdout)


def rescaleSum(bedGraph: str, chromSizes: str, window: int,
               regions: str, distanceTransform: bool, includeZero: bool):
    """ Rescale intervals to a set window size and sum interval scores
        within a window. Optionally performing distancing correction. """

    scores = readBedgraphSum(bedGraph, window)
    template = makeTemplate(
        window, mode, regions, chromSizes,
        distanceTransform, includeZero, threshold=None)

    json.dump(makeRescaledSum(scores, template), sys.stdout)


def makeRescaled(scores, template):
    """ Generate rescaled dictionary for binary and count mode. """

    for chrom, size in template['chromSizes'].items():
        for start in range(0, size, template['window']):
            # Skip windows where start not in regions
            if not validRegion(template['regions'], chrom, start):
                continue
            try:
                score = scores[chrom][start]
            except KeyError:
                score = 0
            if template['mode'] == 'binary':
                score = score > 0.5
            else:
                score = round(score)
            template['data'][chrom][start] = score
    return template


def makeRescaledSum(scores, template):
    """ Generate rescaled dictionary for binary and count mode. """

    for chrom, size in template['chromSizes'].items():
        # Reset prevScore for each chromosome
        prevScore = None
        for start in range(0, size, template['window']):
            # Skip windows where start not in regions
            if not validRegion(regions, chrom, start):
                prevScore = None
                continue
            try:
                score = scores[chrom][start]
            except KeyError:
                prevScore = 0
                if template['includeZero']:
                    score = 0
                else:
                    continue
            if template['distanceTransform']:
                try:
                    template['data'][chrom][start] = score - prevScore
                except TypeError:
                    pass  # Skip windows where prevScore set to None
            else:
                template['data'][chrom][start] = score
            prevScore = score
    return template


def makeTemplate(window, mode, regions, chromSizes,
                 distanceTransform, includeZero, threshold):
    """ Return template rescaled JSON with metadata """
    template = {'window':            window,
                'mode':              mode,
                'threshold':         threshold,
                'regions':           readRegions(regions),
                'chromSizes':        readChromSizes(chromSizes),
                'distanceTransform': distanceTransform,
                'includeZero':       includeZero,
                'data':              defaultdict(dict)}
    return template


def readBedgraph(bedGraph, window, threshold=None):
    scores = defaultdict(dict)
    with open(bedGraph) as fh:
        for line in fh:
            chrom, start, end = splitPos(line)
            regionLength = end - start
            score = score > threshold if threshold is not None else 1
            for base in range(start, end):
                pos = getWindow(base, window)
                try:
                    scores[chrom][pos] += score / regionLength
                except KeyError:
                    scores[chrom][pos] = score / regionLength
    return scores


def readBedgraphSum(bedGraph, window):
    scores = defaultdict(dict)
    with open(bedGraph) as fh:
        for i, line in enumerate(fh):
            if bed:
                chrom, start, end, score = splitBed(line)
            else:  # Input is bedgraph format
                chrom, start, end, score = splitBedgraph(line)
            regionLength = end - start
            for base in range(start, end):
                pos = getWindow(base, window)
                try:
                    scores[chrom][pos] += score / regionLength
                except KeyError:
                    scores[chrom][pos] = score / regionLength
    return scores


def splitPos(line):
    chrom, start, end = line.split()[:3]
    return chrom, int(start), int(end)


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
            regions[chrom].append(range(start, end))
    return regions


def validRegion(regions, chrom, start):
    """ Return True if interval present in regions dict """
    # All regions assumed valid if no regions dict provided
    if regions is None:
        return True
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
    except ValueError as e:
        logging.exception('Input BED file does not contain score column.')
    return chrom, int(start), int(end), float(score)


def getWindow(pos, window):
    """ Return 0-based window start for a given position"""
    return (pos // window) * window


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
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    subparser = parser.add_subparsers(
        title='required commands', description='',
        metavar='Commands', help='Description:')

    baseParser = argparse.ArgumentParser(add_help=False)
    requiredNamed = baseParser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--window', type=int, required=True,
        help='Rescaled interval window size.')
    baseParser.add_argument(
        '--regions', help='BED file indicating regions to process.')
    chromSizeParser = argparse.ArgumentParser(add_help=False)
    chromSizeParser.add_argument(
        'chromSizes', help='Chromosome sizes file.')
    infileParser = argparse.ArgumentParser(add_help=False)
    infileParser.add_argument(
        'bedGraph', help='BedGraph/BED interval file to rescale.')
    infileParser2 = argparse.ArgumentParser(add_help=False)
    infileParser2.add_argument(
        'bedGraph',
        help='BedGraph/BED interval file to rescale. '
             'If BED must also set "--bed".')

    count = subparser.add_parser(
        'count',
        description='Rescale intervals to a set window size and '
                    'count number of intervals falling into each '
                    'window. Interval scores ignored.',
        help='Count number of intervals within each window.',
        epilog=parser.epilog,
        parents=[mainParent, infileParser, chromSizeParser, baseParser])
    count.set_defaults(function=rescaleBinCount)

    binary = subparser.add_parser(
        'binary',
        description='Rescale intervals to a set window size and '
                    'assign each window True/False depending on '
                    'if it contains atleast 1 interval. Interval '
                    'scores ignored.',
        help='Set True/False if window contains atleast 1 interval.',
        epilog=parser.epilog,
        parents=[mainParent, infileParser, chromSizeParser, baseParser])
    binary.add_argument(
        '--threshold', type=float, default=False,
        help='Score threshold for determining binary '
             'intervals (default: None')
    binary.set_defaults(function=rescaleBinCount)

    sum = subparser.add_parser(
        'sum',
        description=rescaleSum.__doc__,
        help='Sum scores across intervals in a window.',
        epilog=parser.epilog,
        parents=[mainParent, infileParser2, chromSizeParser, baseParser])
    sum.add_argument(
        '--includeZero', action='store_true',
        help='Include windows with a 0 score (default: %(default)s)')
    sum.add_argument(
        '--distanceTransform', action='store_true',
        help='Perform differencing to remove series '
             'dependence. (default: %(default)s)')
    sum.add_argument(
        '--bed', action='store_true',
        help='Treat input as BED format (default: %(default)s)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
