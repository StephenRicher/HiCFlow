#!/usr/bin/env python3

""" Summarise bedgraph score per interval region """

import sys
import logging
import argparse
import numpy as np
import pandas as pd
from typing import List
from collections import defaultdict
from utilities import setDefaults, createMainParent
try:
    from pandarallel import pandarallel
except ModuleNotFoundError:
    pass

__version__ = '1.0.0'


def scoreAllIntervals(
        bedGraph: str, beds: List, summary: str,
        groupName: bool, multiplier: int, threads: int):
    intervalSize, shift, bedGraph = readBedgraph(bedGraph)
    bed = pd.concat(readBed(bed) for bed in beds)
    if (threads > 1) and ('pandarallel' in sys.modules):
        pandarallel.initialize(nb_workers=threads, verbose=0)
        bed[['score', 'std']] = bed.parallel_apply(
            scoreInterval, args=(bedGraph, intervalSize, shift), axis=1)
    else:
        bed[['score', 'std']] = bed.apply(
            scoreInterval, args=(bedGraph, intervalSize, shift), axis=1)
    # Normalise score by region length
    bed['score'] = (bed['score'] / (bed['end'] - bed['start'])) * multiplier
    bed.to_csv(sys.stdout, sep='\t', index=False, header=False)
    if not groupName:
        bed['name'] = 'data'
    bed.groupby('name')['score'].describe().to_csv(summary, sep='\t')


def scoreInterval(bed, bedGraph, intervalSize, shift):
    binCount = countBins(
        np.arange(bed['start'], bed['end']), intervalSize, shift)
    score = []
    for pos, count in binCount.items():
        try:
            score.extend([bedGraph[bed['chrom']][pos]] * count)
        except KeyError:
            pass
    score = np.array(score)
    scoreSum = score.sum()
    scoreStd = np.nan if len(score) == 0 else score.std()
    return pd.Series([scoreSum, scoreStd])


def countBins(positions, binSize: int, shift: int):
    """ Return dictionary number of positions mapping to each bin """
    bins = (positions - (positions % binSize)) + shift
    unique, counts = np.unique(bins, return_counts=True)
    return dict(zip(unique, counts))


def readBedgraph(file: str):
    """ Read bedgraph into nested dictionary, keyed
        by start position of bin. """
    bedGraph = defaultdict(dict)
    intervals = set()
    with open(file) as fh:
        for line in fh:
            chrom, start, end, score = line.split()
            intervalSize = int(end) - int(start)
            normScore = float(score) / intervalSize
            # Binsizes don't always align with 0 e.g. 2, 1002, 2002,
            shift = int(start) % intervalSize
            intervals.add(intervalSize)
            if start in bedGraph[chrom]:
                logging.error('Bedgraph contains overlapping intervals.')
                raise ValueError
            if len(intervals) > 1:
                logging.error('Bedgraph is not constant interval size.')
                raise ValueError
            bedGraph[chrom][int(start)] = normScore
    return intervalSize, shift, bedGraph


def addName(bed):
    """ Create BED entry name if absent """
    if pd.isna(bed['name']):
        return f'{bed["chrom"]}:{bed["start"]}-{bed["end"]}'
    else:
        return bed['name']


def readBed(file: str):
    """ Read BED file into Pandas """
    columns = {'chrom': str, 'start': int, 'end': int,
               'name': str, 'score': str, 'strand': str}
    bed = pd.read_csv(
        file, sep='\t', comment='#', header=None,
        dtype=columns, names=columns.keys())
    bed['name'] = bed.apply(addName, axis=1)
    return bed


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument(
        'bedGraph', help='BedGraph (constant interval) with scores.')
    parser.add_argument(
        'beds', metavar='BED', nargs='+',
        help='BED files indicating regions to process.')
    parser.add_argument(
        '--summary', type=str, default=sys.stderr,
        help='File to write score summary statistics  (default: stderr')
    parser.add_argument(
        '--threads', type=int, default=1,
        help='Threads for parallel processing (default: %(default)s)')
    parser.add_argument(
        '--multiplier', type=int, default=1,
        help='Compute score per unit base pair (default: %(default)s)')
    parser.add_argument(
        '--groupName', action='store_true',
        help='Output summary statistics grouped by BED name '
             '(default: %(default)s)')
    parser.set_defaults(function=scoreAllIntervals)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
