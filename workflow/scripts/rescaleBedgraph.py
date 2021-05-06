#!/usr/bin/env python3

""" Convert bedgraph intervals to new binsize and write to PANDAS. """

import sys
import json
import logging
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
from utilities import setDefaults, createMainParent, readChromSizes
from bedgraphUtils import splitScore, readRegions, splitPos


__version__ = '1.0.0'


def rescaleBinCount(bedGraph: str, chromSizes: str, out: str, binSize: int,
                    name: str, regions: str, filetype: str = None,
                    threshold: float = False):

    if (threshold) and (filetype is None):
        logging.error(
            'If --threshold is set then --filetype must also be specified')
        return 1
    mode = 'count'
    scores = readBedgraph(bedGraph, binSize, mode, filetype, threshold)
    name = bedGraph if name is None else name
    template = makeTemplate(
        name, binSize, mode, regions, chromSizes,
        distanceTransform=False, includeZero=True, threshold=threshold)
    df = convert2pandas(makeRescaledSum(scores, template))
    df.to_pickle(out)


def rescaleSum(bedGraph: str, chromSizes: str, out: str, binSize: int, name: str,
               regions: str, filetype: str, distanceTransform: bool,
               includeZero: bool):
    """ Rescale intervals to a set binsize and sum interval scores
        within a bin. Optionally performing distancing correction. """

    mode = 'sum'
    scores = readBedgraph(bedGraph, binSize, mode, filetype)
    name = bedGraph if name is None else name
    template = makeTemplate(
        name, binSize, mode, regions, chromSizes,
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



def makeRescaledSum(scores, template):
    """ Generate rescaled dictionary for binary and count mode. """

    for chrom, size in template['meta']['chromSizes'].items():
        # Reset prevScore for each chromosome
        prevScore = None
        for start in range(0, size, template['meta']['binSize']):
            # Skip bins where start not in regions
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
                    pass  # Skip bins where prevScore set to None
            else:
                template['data'][chrom][start] = score
            prevScore = score
    return template


def makeTemplate(name, binSize, mode, regions, chromSizes,
                 distanceTransform, includeZero, threshold):
    """ Return template rescaled JSON with metadata """
    template = {'data':              defaultdict(dict),
                'meta':
                    {'name':              name                      ,
                     'binSize':           binSize                  ,
                     'mode':              mode                      ,
                     'threshold':         threshold                 ,
                     'regions':           readRegions(regions)      ,
                     'chromSizes':        readChromSizes(chromSizes),
                     'distanceTransform': distanceTransform         ,
                     'includeZero':       includeZero               ,}
                }
    return template


def readBedgraph(bedGraph, binSize, mode, filetype=None, threshold=None):
    scores = defaultdict(dict)
    assert mode in ['sum', 'count']
    with open(bedGraph) as fh:
        for i, line in enumerate(fh):
            if mode == 'sum':
                chrom, start, end, score = splitScore(line, filetype)
            else:
                if isinstance(threshold, float):
                    chrom, start, end, score = splitScore(line, filetype)
                    score = score > threshold
                else:
                    chrom, start, end = splitPos(line)
                    score = 1
            if score == 0:
                continue
            score = score / (end - start)
            perBaseBin = getBin(np.array(range(start, end)), binSize)
            bin, counts = np.unique(perBaseBin, return_counts=True)
            counts = counts * score
            perBinScore = dict(zip(bin, counts))
            for bin, binScore in perBinScore.items():
                try:
                    scores[chrom][bin] += binScore
                except KeyError:
                    scores[chrom][bin] = binScore
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


def getBin(pos, binSize):
    """ Return 0-based bin start for a given position"""
    return (pos // binSize) * binSize


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
        '--regions',
        help='BED file indicating regions to process.')

    formatParser = argparse.ArgumentParser(add_help=False)
    formatRequired = formatParser.add_argument_group('required named arguments')
    formatRequired.add_argument(
        '--filetype',  required=True, choices=['bed', 'bedgraph'],
        help='Input filetype to correctly retrieve score column.')

    chromSizeParser = argparse.ArgumentParser(add_help=False)
    chromSizeParser.add_argument(
        'chromSizes', help='Chromosome sizes file.')

    infileParser = argparse.ArgumentParser(add_help=False)
    infileParser.add_argument(
        'bedGraph', help='BedGraph/BED interval file to rescale.')

    binSizeParser = argparse.ArgumentParser(add_help=False)
    binSizeRequired = binSizeParser.add_argument_group('required named arguments')
    binSizeRequired.add_argument(
        '--binSize', type=int, required=True,
        help='Rescaled interval binSize size.')

    count = subparser.add_parser(
        'count',
        description='Rescale intervals to a set binsize and '
                    'count number of intervals falling into each '
                    'bin. Interval scores ignored.',
        help='Count number rescaleSumof intervals within each bin.',
        epilog=parser.epilog,
        parents=[mainParent, infileParser, chromSizeParser,
                 baseParser, regionsParser, binSizeParser])
    count.add_argument(
        '--threshold', type=float, default=False,
        help='Minimum score threshold for determining binary '
             'intervals (default: None)')
    count.add_argument(
        '--filetype',  choices=['bed', 'bedgraph'],
        help='Input format to correctly retrieve score column. '
             'Only required if --threshold used (default: %(default)s)')
    count.set_defaults(function=rescaleBinCount)

    sum = subparser.add_parser(
        'sum',
        description=rescaleSum.__doc__,
        help='Sum scores across intervals in a bin.',
        epilog=parser.epilog,
        parents=[mainParent, infileParser, chromSizeParser,
                 baseParser, regionsParser, binSizeParser, formatParser])
    sum.add_argument(
        '--includeZero', action='store_true',
        help='Include bins with a 0 score (default: %(default)s)')
    sum.add_argument(
        '--distanceTransform', action='store_true',
        help='Perform differencing to remove series '
             'dependence. (default: %(default)s)')
    sum.set_defaults(function=rescaleSum)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
