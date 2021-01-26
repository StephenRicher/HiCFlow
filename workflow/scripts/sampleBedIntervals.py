#!/usr/bin/env python3

""" Randomly create BED intervals from referenceBED based on
    length in sampleBED  """

import os
import sys
import random
import logging
import argparse
from itertools import repeat
from bedgraphUtils import readBedLength
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def sampleIntervals(referenceBed: str, sampleBed: str, sampleFormat: str,
                    nReps: float, seed: float):
    random.seed(seed)
    # Set number of attempts to find suitable position
    maxAttempts = 1000
    regions = readBedLength(referenceBed)
    print(f'# Random intervals sampled from {os.path.basename(referenceBed)} '
          f'matching interval lengths in {os.path.basename(sampleBed)}. '
          f'Total repeat interval sets: {nReps}.')
    for rep in range(nReps):
        nAttempts = 0
        sampleIntervals = readBedLength(sampleBed, fileType=sampleFormat)
        while (len(sampleIntervals) > 0):
            if nAttempts > maxAttempts:
                logging.error(
                    f'Intervals of the following length could not be found '
                    f'within boundaries of the reference:\n {sampleIntervals}')
                break
            nSamples = 0
            repeatIntervals = {}
            # Select BEDS, weight by length, with replacement
            selections = random.choices(
                list(regions.keys()),
                weights=list(regions.values()), k=len(sampleIntervals))
            for selection, interval in zip(selections, sampleIntervals):
                pos = random.choice(selection.interval)
                end = pos + interval.regionLength
                # Interval extends beyond boundary - must repeat
                if end > selection.end:
                    repeatIntervals[interval] = interval.regionLength
                elif sampleFormat == 'links':
                    print(selection.chrom, pos, pos + interval.bin1Length,
                          selection.chrom, end - interval.bin2Length, end, rep,
                          sep='\t')
                else:
                    print(selection.chrom, pos, end, rep, sep='\t')
            nAttempts += 1
            sampleIntervals = repeatIntervals



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
        '--sampleFormat', choices=['bed', 'links'], default='bed',
        help='Input format of sampleBED to determine whether to produce single '
             'or paired intervals.')
    parser.add_argument(
        '--nReps', type=int, default=1,
        help='Number of sets of intervals to generate.')
    parser.add_argument(
        '--seed', default=None, type=float,
        help='Seed for random number generation (default: %(default)s)')
    parser.set_defaults(function=sampleIntervals)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
