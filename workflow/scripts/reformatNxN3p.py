#!/usr/bin/env python3

""" Create new FastQC zip directory with a new sample name for multiQC."""

import os
import sys
import logging
import argparse
import fileinput


__version__ = '1.0.0'


def main(file, binSize, refName, **kwargs):

    with fileinput.input(file) as fh:
        for i, line in enumerate(fh):
            # Skip header line
            if i == 0: continue
            counts = line.split()
            start = (i - 1) * binSize
            end = start + binSize
            # Remove coordinate columns
            counts = counts[2:]
            # Add new columns to start
            counts = [refName, str(start), str(end)] + counts
            print('\t'.join(counts))


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'binSize', type=int,
        help='Bin size of HiC matrix')
    custom.add_argument(
        'refName',
        help='Reference name e.g. chromosome or capture region name.')
    custom.add_argument(
        'file', nargs='?', default=[],
        help='HiC matrix in homer format (default: stdin)')
    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'

    base = argparse.ArgumentParser(add_help=False)
    base.add_argument(
        '--version', action='version', version=f'%(prog)s {__version__}')
    base.add_argument(
        '--verbose', action='store_const', const=logging.DEBUG,
        default=logging.INFO, help='verbose logging for debugging')

    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[base, custom])
    args = parser.parse_args()

    log_format='%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
    logging.basicConfig(level=args.verbose, format=log_format)

    return args


if __name__ == '__main__':
    args = parse_arguments()
    return_code = args.function(**vars(args))
    logging.shutdown()
    sys.exit(return_code)
