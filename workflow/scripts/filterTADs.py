#!/usr/bin/env python3

""" Filter TAD domains file to return only region overlapping defined region """

import os
import sys
import logging
import argparse
import fileinput
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def filterTADs(TADdomains: str, chromRef: str, startRef: int, endRef: int):
    with fileinput.input(TADdomains) as fh:
        for line in fh:
            chrom, start, end = line.split()[:3]
            start = int(start)
            end = int(end)
            # Exclude entries from different chromosome
            if chrom != chromRef:
                continue
            # Exclude entries with start lower than startRef
            if (startRef is not None) and (start < startRef):
                continue
            # Exclude entries with end higher than endRef
            if (endRef is not None) and (end > endRef):
                continue
            print(line, end='')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'TADdomains', nargs='?', default=[],
        help='TAD domains file containing intervals '
             'to be filtered (default: stdin)')
    parser.add_argument(
        '--start', type=int, dest='startRef',
        help='Start position (0-based) of region to '
             'filter (default: %(default)s)')
    parser.add_argument(
        '--end', type=int, dest='endRef',
        help='End position (1-based of region to filter) '
             '(default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--chrom', dest='chromRef',
        help='Chromosome of filtered region (default: %(default)s)')
    parser.set_defaults(function=filterTADs)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
