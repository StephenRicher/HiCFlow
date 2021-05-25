#!/usr/bin/env python3

""" Convert bedgraph intervals to new binsize and write to PANDAS. """

import sys
import json
import logging
import argparse
import itertools
import numpy as np
import pandas as pd
from collections import defaultdict
from utilities import setDefaults, createMainParent, readChromSizes


__version__ = '1.0.0'


def rescaleCount(bedGraph: str, chromSizes: str, out: str, binSize: int,
                 name: str, regions: str, minOverlap: int,
                 filetype: str = None, threshold: float = False):
    chromSizes = readChromSizes(chromSizes)
    if threshold:
        if filetype is None:
            logging.error('If --threshold is set then --filetype '
                          'must also be specified')
            return 1
        # Select coordinate and score columns depending on filetype
        bed = readBedgraph(bedGraph, filetype=filetype, score=True)
        # Filter score regions by threshold
        bed = bed.loc[bed['score'] >= threshold]
    else:
        bed = readBedgraph(bedGraph, score=False)
    rescaledByChrom = []
    for chrom, df in bed.groupby('chrom'):
        # For count mode, region scores are set to the region length
        df['score'] = df['end'] - df['start']
        # Split regions into seperate rows according to their bin
        df = rescaleInterval(df, binSize)
        # Count a region if atleast X number of bases overlap a bin
        df['count'] = df['overlap'] > minOverlap
        # Start position of last bin in chromosome
        maxBin = (chromSizes[chrom] // binSize) * binSize
        df = (df
            .groupby('startBin').apply(baseOverlap, binSize)
            .rename({0: 'overlap', 1: 'count'}, axis=1)
            .reindex(range(0, maxBin, binSize), fill_value=0)
            .reset_index()
            .rename({'startBin': 'start'}, axis=1))
        df['chrom'] = chrom
        rescaledByChrom.append(df)
    # Merge across chromosomes and save
    try:
        rescaledByChrom = pd.concat(rescaledByChrom, axis=0)
    except ValueError:
        names = {'start': int, 'overlap': int, 'count': int, 'chrom': str}
        rescaledByChrom = pd.DataFrame(columns=names.keys()).astype(names)
    rescaledByChrom['end'] = rescaledByChrom['start'] + binSize
    rescaledByChrom.attrs['name'] = name
    rescaledByChrom.attrs['mode'] = 'count'
    cols = ['chrom', 'start', 'end', 'overlap', 'count']
    rescaledByChrom[cols].to_pickle(out)



def rescaleSum(bedGraph: str, chromSizes: str, out: str, binSize: int,
               name: str, regions: str, filetype: str):
    """ Rescale intervals to a set binsize and sum interval scores
        within a bin. Optionally performing distancing correction. """
    chromSizes = readChromSizes(chromSizes)
    bed = readBedgraph(bedGraph, filetype=filetype, score=True)
    rescaledByChrom = []
    for chrom, df in bed.groupby('chrom'):
        # Split regions into seperate rows according to their bin
        df = rescaleInterval(df, binSize)
        # Start position of last bin in chromosome
        maxBin = (chromSizes[chrom] // binSize) * binSize
        df = (df
            .groupby('startBin')['score']
            .sum()
            .rename('score')
            .reindex(range(0, maxBin, binSize), fill_value=0)
            .reset_index()
            .rename({'startBin': 'start'}, axis=1))
        df['chrom'] = chrom
        df['scoreDiff'] = df['score'].diff(1)
        rescaledByChrom.append(df)
    # Merge across chromosomes and save
    rescaledByChrom = pd.concat(rescaledByChrom, axis=0)
    rescaledByChrom['end'] = rescaledByChrom['start'] + binSize
    cols = ['chrom', 'start', 'end', 'score', 'scoreDiff']
    rescaledByChrom.attrs['name'] = name
    rescaledByChrom.attrs['mode'] = 'sum'
    rescaledByChrom[cols].to_pickle(out)


def rescaleInterval(bed, binSize):
    df2 = pd.DataFrame(
        [x for x in itertools.chain.from_iterable(
            [IntervGen(row.start, row.end, row.score, binSize) for row in bed.itertuples()])],
        columns=['start', 'end', 'score'])
    # Calculate number of bases each region overlaps a particular bin
    df2['overlap'] = df2['end'] - df2['start']
    # Round start positions down to nearest multiple of binSize
    df2['startBin'] = (df2['start'] // binSize) * binSize
    return df2


def readBedgraph(file, filetype='bedgraph', score=True):
    if score == False:
        useCols = [0, 1, 2]
        dtypes = {'chrom': str, 'start': int, 'end': int}
    else:
        useCols = [0, 1, 2, 4] if filetype == 'bed' else [0, 1, 2, 3]
        dtypes = {'chrom': str, 'start': int, 'end': int, 'score': float}
    return pd.read_csv(
        file, usecols=useCols, comment='#', header=None,
        names=dtypes.keys(), dtype=dtypes, sep='\t')


def IntervGen(st, en, val, binSize):
    """ https://stackoverflow.com/questions/62651415/change-interval-length-in-pandas """
    nxtSt = (st // binSize + 1) * binSize
    totalSize = en - st
    origVal = val
    while st < en:
        if en <= nxtSt:
            yield st, en, val
            return
        currSize = nxtSt - st
        currVal = origVal * currSize / totalSize
        yield st, nxtSt, currVal
        st = nxtSt
        nxtSt += binSize
        val -= currVal


def baseOverlap(df, binSize):
    """ Proportion of non-overlapping bases in a bin. """
    allPos = []
    for start, end, startBin in zip(df['start'], df['end'], df['startBin']):
        end = min(end, startBin + binSize)
        allPos.extend(list(range(start, end)))
    count = df['count'].sum()
    return pd.Series([len(set(allPos)) / binSize, count])


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
        help='Count number of intervals within each bin.',
        epilog=parser.epilog,
        parents=[mainParent, baseParser, regionsParser,
                 binSizeParser, chromSizeParser, infileParser])
    count.add_argument(
        '--threshold', type=float, default=False,
        help='Minimum score threshold for determining binary '
             'intervals (default: None)')
    count.add_argument(
        '--minOverlap', default=1000, type=int,
        help='Minimum number of bases required to overlap bin to be counted '
             'as present in that bin (default: %(default)s)')
    count.add_argument(
        '--filetype',  choices=['bed', 'bedgraph'],
        help='Input format to correctly retrieve score column. '
             'Only required if --threshold used (default: %(default)s)')
    count.set_defaults(function=rescaleCount)

    sum = subparser.add_parser(
        'sum',
        description=rescaleSum.__doc__,
        help='Sum scores across intervals in a bin.',
        epilog=parser.epilog,
        parents=[mainParent, baseParser, regionsParser, binSizeParser,
                 formatParser, chromSizeParser, infileParser])
    sum.set_defaults(function=rescaleSum)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
