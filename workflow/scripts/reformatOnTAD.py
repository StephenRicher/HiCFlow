#!/usr/bin/env python3

""" Reformat OnTAD output to BED coordinates."""

import sys
import argparse
import fileinput
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def reformatOnTAD(ontad: str, scale: int, chrom: str, binSize: int, maxPos: int):
    assert scale >= 0
    with fileinput.input(ontad) as fh:
        for i, line in enumerate(fh):
            # Skip first line
            if i == 0:
                continue
            # Start and End bins are 1-based in OnTAD
            start, end = line.split()[:2]
            score = line.split()[4]
            # Retrieve midpoint coordinates of each interval
            startMid = int((int(start) * binSize) - (binSize / 2))
            # Rest end position to maximum position.
            # Usually if chromosome end is within a bin.
            endPos = min(int(end) * binSize, maxPos)
            endMid = int(endPos - (binSize / 2))
            # Shift coordinates as required
            startMid = str(startMid  + scale)
            endMid = str(endMid + scale)
            print(chrom, startMid, endMid, '.', score, '.', startMid, endMid, '211,211,211', sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])

    parser.add_argument(
        'ontad', nargs='?', default=[],
        help='Result from OnTAD (default: stdin)')
    parser.add_argument(
        '--scale', type=int, default=0,
        help='0-based genomic coordinate of start position '
             'to rescale (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--chrom', required=True, help='Chromosome of OnTAD result.')
    requiredNamed.add_argument(
        '--binSize', required=True, type=int,
        help='Binsize used for OnTAD calculation.')
    requiredNamed.add_argument(
        '--maxPos', required=True, type=int,
        help='Max position of matrix.')

    parser.set_defaults(function=reformatOnTAD)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
