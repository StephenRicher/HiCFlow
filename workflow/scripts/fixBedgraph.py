#!/usr/bin/env python3

""" Fix final position of bedgraph entry """


import sys
import argparse
import fileinput
from utilities import setDefaults

__version__ = '1.0.0'


def fixBedgraph(file: str, pos: int):

    with fileinput.input(file) as fh:
        while True:
            try:
                line = next(fh).strip()
                print(prevLine)
            except StopIteration:
                # Last line in file
                ref, start, end, score = prevLine.split()
                print(ref, start, pos, score, sep='\t')
                break
            except UnboundLocalError:
                pass
            prevLine = line


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'file', nargs='?',
        help='Input bedgraph.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--pos', type=int, required=True,
        help='Value to assign to End position of last bedgraph entry.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(fixBedgraph(**vars(args)))
