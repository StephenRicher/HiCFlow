#!/usr/bin/env python3

""" Parse HapCUT2 block file to extract block with most phased sites. """

import os
import sys
import logging
import argparse


__version__ = '1.0.0'


def main(file, **kwargs):

    try:
        bestBlockLine = getBestBlockLine(file)
    except UnboundLocalError:
        logging.warning(f'{file} is empty.')
        return 0
        
    with open(file) as fh:
        for i, line in enumerate(fh):
            if i == bestBlockLine:
                while True:
                    try:
                        line = next(fh)
                    except StopIteration:
                        break
                    if line.startswith('********'):
                        break
                    ref, pos = line.split()[3:5]
                    print(ref, pos, pos, sep='\t')
                break


def getBestBlockLine(file):
    """ Return line number of start of best block """

    bestPhased = 0
    with open(file) as f:
        for i, line in enumerate(f):
            if line.startswith('BLOCK'):
                nPhased = int(line.split()[6])
                if nPhased > bestPhased:
                    bestBlockLine = i
                    bestPhased = nPhased
    return bestBlockLine


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument('file', help='HapCUT2 block file. )')
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
