#!/usr/bin/env python3

""" Reformat OnTAD output BED file to links scaled to correct coordinates."""

import sys
import logging
import argparse
import fileinput


__version__ = '1.0.0'


def main(start, bin, file, **kwargs):

    with fileinput.input(file) as fh:
        for i, line in enumerate(fh):
            # Skip header line
            if i == 0:
                continue
            columns = line.split()
            chr = columns[0]
            # Lower bound of interaction
            lower = int(columns[1]) + start
            # Upper bound of interaction
            upper = int(columns[2]) + start

            print(chr, lower, lower+bin, chr, upper-bin, upper, 0, sep='\t')


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'start', type=int,
        help='Genomic coordinates of start position (0 based).')
    custom.add_argument(
        'bin', type=int,
        help='Bin size of HiC matrix.')
    custom.add_argument(
        'file', nargs='?', default=[],
        help='OnTAD output in BED format (default: stdin)')
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
