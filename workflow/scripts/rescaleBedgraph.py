#!/usr/bin/env python3

""" Convert bedgraph intervals to new window size and write to JSON. """

import sys
import json
import logging
import argparse
import pandas as pd
from collections import defaultdict
from utilities import setDefaults, createMainParent
from bedgraphUtils import splitScore, readRegions, splitPos


__version__ = '1.0.0'


def rescaleBinCount(bedGraph: str, chromSizes: str, out: str, window: int,
                    name: str, regions: str, format: str = None,
                    threshold: float = None):

    if threshold is not None:
        mode = 'binary'
        if (threshold is not False) and (format is None):
            logging.error(
                'If --threshold is set then --format must also be specified')
            return 1
    else:
        mode = 'count'

    scores = readBedgraph(bedGraph, window, format, threshold)
    name = bedGraph if name is None else name
    template = makeTemplate(
        name, window, mode, regions, chromSizes,
        distanceTransform=False, includeZero=True, threshold=threshold)
    df = convert2pandas(makeRescaledSum(scores, template))
    df.to_pickle(out)


def rescaleSum(bedGraph: str, chromSizes: str, out: str, window: int, name: str,
               regions: str, format: str, distanceTransform: bool,
               includeZero: bool):
    """ Rescale intervals to a set window size and sum interval scores
        within a window. Optionally performing distancing correction. """

    scores = readBedgraphSum(bedGraph, window, format)
    name = bedGraph if name is None else name
    mode = 'sum'
    template = makeTemplate(
        name, window, mode, regions, chromSizes,
        distanceTransform, includeZero, threshold=None)
    df = convert2pandas(makeRescaledSum(scores, template))
    df.to_pickle(out)


def convert2pandas(template):
    """ Write dictionary template to pandas DF with metadata """
    columnIndex = pd.MultiIndex.from_tuples(
        [(template['meta']['mode'], template['meta']['name'])],
        names=('mode', 'name'))
    df = pd.DataFrame.from_dict(
        template['data'], orient='index').stack().to_frame()
    df.columns = columnIndex
    df = df.reindex(df.index.set_names(['chromosome', 'start']))
    df.attrs['meta'] = template['meta']
    return df


def makeRescaled(scores, template):
    """ Generate rescaled dictionary for binary and count mode. """

    for chrom, size in template['meta']['chromSizes'].items():
        for start in range(0, size, template['meta']['window']):
            # Skip windows where start not in regions
            if not validRegion(template['meta']['regions'], chrom, start):
                continue
            try:
                score = scores[chrom][start]
            except KeyError:
                score = 0
            if template['meta']['mode'] == 'binary':
                score = score > 0.5
            else:
                score = round(score)
            template['data'][chrom][start] = score
    return template


def makeRescaledSum(scores, template):
    """ Generate rescaled dictionary for binary and count mode. """

    for chrom, size in template['meta']['chromSizes'].items():
        # Reset prevScore for each chromosome
        prevScore = None
        for start in range(0, size, template['meta']['window']):
            # Skip windows where start not in regions
            if not validRegion(template['meta']['regions'], chrom, start):
                prevScore = None
                continue
            try:
                score = scores[chrom][start]
            except KeyError:
                prevScore = 0
                if template['meta']['includeZero']:
                    score = 0
                else:
                    continue
            if template['meta']['distanceTransform']:
                try:
                    template['data'][chrom][start] = score - prevScore
                except TypeError:
                    pass  # Skip windows where prevScore set to None
            else:
                template['data'][chrom][start] = score
            prevScore = score
    return template


def makeTemplate(name, window, mode, regions, chromSizes,
                 distanceTransform, includeZero, threshold):
    """ Return template rescaled JSON with metadata """
    template = {'data':              defaultdict(dict),
                'meta':
                    {'name':              name                      ,
                     'window':            window                    ,
                     'mode':              mode                      ,
                     'threshold':         threshold                 ,
                     'regions':           readRegions(regions)      ,
                     'chromSizes':        readChromSizes(chromSizes),
                     'distanceTransform': distanceTransform         ,
                     'includeZero':       includeZero               ,}
                }
    return template


def readBedgraph(bedGraph, window, format=None, threshold=None):
    scores = defaultdict(dict)
    with open(bedGraph) as fh:
        for line in fh:
            # Score column must be retrieved if threshold set
            if isinstance(threshold, float):
                chrom, start, end, score = splitScore(line, format)
                score > threshold
            else:
                chrom, start, end = splitPos(line)
                score = 1
            regionLength = end - start
            for base in range(start, end):
                pos = getWindow(base, window)
                try:
                    scores[chrom][pos] += score / regionLength
                except KeyError:
                    scores[chrom][pos] = score / regionLength
    return scores


def readBedgraphSum(bedGraph, window, format):
    scores = defaultdict(dict)
    with open(bedGraph) as fh:
        for i, line in enumerate(fh):
            chrom, start, end, score = splitScore(line, format)
            regionLength = end - start
            for base in range(start, end):
                pos = getWindow(base, window)
                try:
                    scores[chrom][pos] += score / regionLength
                except KeyError:
                    scores[chrom][pos] = score / regionLength
    return scores


def validRegion(regions, chrom, start):
    """ Return True if interval present in regions dict """
    # All regions assumed valid if no regions dict provided
    if regions is None:
        return True
    for interval in regions[chrom]:
        if start in interval:
            return True
    return False


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
        '--out', required=True,
        help='Path to save pickled dataframe.')
    baseParser.add_argument(
        '--name', help='Name to store in metadata. Defaults to infile path.')

    regionsParser = argparse.ArgumentParser(add_help=False)
    regionsRequired = regionsParser.add_argument_group('required named arguments')
    regionsRequired.add_argument(
        '--regions', required=True,
        help='BED file indicating regions to process.')

    formatParser = argparse.ArgumentParser(add_help=False)
    formatRequired = formatParser.add_argument_group('required named arguments')
    formatRequired.add_argument(
        '--format',  required=True, choices=['bed', 'bedgraph'],
        help='Input format to correctly retrive score column.')

    chromSizeParser = argparse.ArgumentParser(add_help=False)
    chromSizeParser.add_argument(
        'chromSizes', help='Chromosome sizes file.')

    infileParser = argparse.ArgumentParser(add_help=False)
    infileParser.add_argument(
        'bedGraph', help='BedGraph/BED interval file to rescale.')

    windowParser = argparse.ArgumentParser(add_help=False)
    windowRequired = windowParser.add_argument_group('required named arguments')
    windowRequired.add_argument(
        '--window', type=int, required=True,
        help='Rescaled interval window size.')

    count = subparser.add_parser(
        'count',
        description='Rescale intervals to a set window size and '
                    'count number of intervals falling into each '
                    'window. Interval scores ignored.',
        help='Count number of intervals within each window.',
        epilog=parser.epilog,
        parents=[mainParent, infileParser, chromSizeParser,
                 baseParser, regionsParser, windowParser])
    count.set_defaults(function=rescaleBinCount)

    binary = subparser.add_parser(
        'binary',
        description='Rescale intervals to a set window size and '
                    'assign each window True/False depending on '
                    'if it contains atleast 1 interval. Interval '
                    'scores ignored.',
        help='Set True/False if window contains atleast 1 interval.',
        epilog=parser.epilog,
        parents=[mainParent, infileParser, chromSizeParser,
                 baseParser, regionsParser, windowParser])
    binary.add_argument(
        '--threshold', type=float, default=False,
        help='Minimum score threshold for determining binary '
             'intervals (default: None)')
    binary.add_argument(
        '--format',  choices=['bed', 'bedgraph'],
        help='Input format to correctly retrieve score column. '
             'Only required if --threshold used (default: %(default)s)')
    binary.set_defaults(function=rescaleBinCount)

    sum = subparser.add_parser(
        'sum',
        description=rescaleSum.__doc__,
        help='Sum scores across intervals in a window.',
        epilog=parser.epilog,
        parents=[mainParent, infileParser, chromSizeParser,
                 baseParser, regionsParser, windowParser, formatParser])
    sum.add_argument(
        '--includeZero', action='store_true',
        help='Include windows with a 0 score (default: %(default)s)')
    sum.add_argument(
        '--distanceTransform', action='store_true',
        help='Perform differencing to remove series '
             'dependence. (default: %(default)s)')
    sum.set_defaults(function=rescaleSum)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
