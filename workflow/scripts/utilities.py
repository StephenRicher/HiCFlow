#!/usr/bin/env python3

import sys
import logging
import argparse

def setDefaults(parser):
    """ Set logging config and return args with associated function """

    args = parser.parse_args()
    logFormat='%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
    try:
        logging.basicConfig(level=args.verbose, format=logFormat)
        del args.verbose
    except AttributeError:
        logging.basicConfig(format=logFormat)
        pass
    try:
        function = args.function
        del args.function
    except AttributeError:
        parser.print_help()
        sys.exit()

    return args, function


def createMainParent(verbose=True, version=None):
    """ Create parser of verbose/version to be added to parser/subparsers. """
    parent = argparse.ArgumentParser(add_help=False)
    if version:
        parent.add_argument('--version', action='version',
            version='%(prog)s {}'.format(version))
    if verbose:
        parent.add_argument(
            '--verbose', action='store_const', const=logging.DEBUG,
            default=logging.INFO, help='verbose logging for debugging')
    return parent


def readChromSizes(file):
    """ Read chromosome sizes to dict """
    chromSizes = {}
    with open(file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chrom, size = line.split()
            chromSizes[str(chrom)] = int(size)
    return chromSizes
