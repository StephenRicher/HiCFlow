#!/usr/bin/env python3

""" Randomly create links intervals within BED input file """

import sys
import random
import logging
import argparse
from bedgraphUtils import Bed
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def sampleLoops(bed: str, links: str, nIntervals: int, length: int, seed: float):
    random.seed(seed)
    allEntries = readBed(bed)
    linkLengths = readLinks(links)
    for linkLength in linkLengths:
        # Create copy of nIntervals value as this is modified per loop
        nIntervals_tmp = nIntervals
        while nIntervals_tmp > 0:
            # Select BEDS, weight by length, with replacement
            selections = random.choices(
                list(allEntries.keys()),
                weights=list(allEntries.values()), k=nIntervals_tmp)
            # Reset count to 0 and add any invalid as we loop through selections
            nIntervals_tmp = 0
            for selection in selections:
                start1 = random.choice(selection.interval)
                end1 = start1 + length
                start2 = start1 + linkLength
                end2 = start2 + length
                if end2 < selection.end:
                    print(selection.chrom, start1, end1, selection.chrom, start2, end2, sep='\t')
                else:
                    # If loop extend beyond the length of the chromosome need to reroll
                    nIntervals_tmp += 1


def readBed(file):
    """ Read BED as dictory of BED (keys) and BED length (values) """
    allEntries = {}
    with open(file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            entry = Bed(line)
            allEntries[entry] = entry.regionLength
    return allEntries


def readLinks(file):
    """ Read links file format and extract unique set of lengths """
    lengths = set()
    with open(file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chr1, start1, end1, chr2, start2, end2 = line.split()[:6]
            length = abs(int(start2) - int(start1))
            lengths.add(length)
    return lengths


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument(
        'bed', metavar='BED',
        help='BED intervals within which regions will be sampled .')
    parser.add_argument(
        'links', metavar='LINKS',
        help='File containing links intervals to determine link lengths.')
    parser.add_argument(
        '--length', default=1, type=int,
        help='The length of the intervals to generate.')
    parser.add_argument(
        '--nIntervals', default=10_000, type=int,
        help='Number of intervals per link distance (default: %(default)s).')
    parser.add_argument(
        '--seed', default=None, type=float,
        help='Seed for random number generation (default: %(default)s)')
    parser.set_defaults(function=sampleLoops)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
