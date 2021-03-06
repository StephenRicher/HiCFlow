#!/usr/bin/env python3

""" Fix final position of bedgraph entry """


import sys
import argparse
import fileinput
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def fixBedgraph(file: str, pos: int):

    with fileinput.input(file) as fh:
        while True:
            try:
                line = next(fh).strip()
                print(prevLine)
            except StopIteration:
                try:
                    # Last line in file
                    ref, start, end, score = prevLine.split()
                    if int(start) < pos:
                        print(ref, start, pos, score, sep='\t')
                except UnboundLocalError:
                    pass # Empty input file
                break
            except UnboundLocalError:
                pass
            prevLine = line


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=fixBedgraph)
    parser.add_argument(
        'file', nargs='?',
        help='Input bedgraph.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--pos', type=int, required=True,
        help='Value to assign to End position of last bedgraph entry.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
