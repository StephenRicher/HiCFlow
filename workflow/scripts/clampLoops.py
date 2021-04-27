#!/usr/bin/env python3

""" Bound bedgraph positions to within a min and max range """


import sys
import argparse
import fileinput
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def clampBedgraph(minPos: int, maxPos: int, file: str):

    with fileinput.input(file) as fh:
        for line in fh:
            ref1, start1, end1, ref2, start2, end2, score = line.strip().split()
            start1 = clamp(int(start1), minPos, maxPos)
            end1 = clamp(int(end1), minPos, maxPos)
            start2 = clamp(int(start2), minPos, maxPos)
            end2 = clamp(int(end2), minPos, maxPos)
            if (end1 > start1) and (end2 > start2):
                print(ref1, start1, end1, ref2, start2, end2, score, sep='\t')


def clamp(x, minimum, maximum):
    return max(minimum, min(x, maximum))


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=clampBedgraph)
    parser.add_argument(
        'minPos', type=int, help='Minimum allowed position in bedgraph.')
    parser.add_argument(
        'maxPos', type=int, help='Maximum allowed position in bedgraph.')
    parser.add_argument('file', nargs='?', help='Input bedgraph.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
