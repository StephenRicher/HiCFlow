#!/usr/bin/env python3

""" Modify cutadpt log file to change sample name read by multiQC """

import sys
import logging
import argparse
import fileinput

__version__ = '1.0.0'


def main(sample, file, **kwargs):

    with fileinput.input(file) as fh:
        for line in fh:
            line = line.strip('\n')
            if 'Command line parameters:' in line:
                line += f' # {sample}'
            print(line)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'sample',
        help='Sample name to append to QC log.')
    custom.add_argument(
        'file', nargs='?', default=[],
        help='Cutadapt QC log (default: stdin)')
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
