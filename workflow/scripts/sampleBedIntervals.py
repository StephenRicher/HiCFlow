#!/usr/bin/env python3

""" Randomly create BED intervals from referenceBED based on
    length in sampleBED  """

import sys
import random
import logging
import argparse
from itertools import repeat
from bedgraphUtils import readBedLength
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def sampleIntervals(referenceBed: str, sampleBed: str, name: str, seed: float):
    random.seed(seed)
    # Set number of attempts to find suitable position
    maxAttempts = 1000
    nAttempts = 0
    regions = readBedLength(referenceBed)
    intervalLengths = list(readBedLength(sampleBed).values())
    name = sampleBed if name is None else name
    while (len(intervalLengths) > 0):
        if nAttempts > maxAttempts:
            logging.error(
                f'Intervals of the following length could not be found within '
                f'the boundaries of the reference: {intervalLengths}')
            break
        nSamples = 0
        repeatIntervals = []
        # Select BEDS, weight by length, with replacement
        selections = random.choices(
            list(regions.keys()),
            weights=list(regions.values()), k=len(intervalLengths))
        for i, selection in enumerate(selections):
            pos = random.choice(selection.interval)
            end = pos + intervalLengths[i]
            # Interval extends beyond boundary - must repeat
            if end > selection.end:
                repeatIntervals.append(intervalLengths[i])
            else:
                print(selection.chrom, pos, end, name, sep='\t')
        nAttempts += 1
        intervalLengths = repeatIntervals



def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument(
        'referenceBed',
        help='BED file containing intervals within '
             'which regions will be sampled .')
    parser.add_argument(
        'sampleBed',
        help='BED file to extract interval lengths of sample.')
    parser.add_argument(
        '--name',
        help='Name to append to BED entry - default to BED filename')
    parser.add_argument(
        '--seed', default=None, type=float,
        help='Seed for random number generation (default: %(default)s)')
    parser.set_defaults(function=sampleIntervals)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
