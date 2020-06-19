#!/usr/bin/env python3

""" Filter interactions of HiCcompare by absolute log fold change and p value.
    Optionally write UP and DOWN interactions to different files.
"""

import sys
import logging
import argparse
import fileinput
from contextlib import ExitStack

__version__ = '1.0.0'


def main(file, up, down, log_fc=0, p_value=0.05, **kwargs):

    with ExitStack() as stack:

        fh = stack.enter_context(fileinput.input(file))
        out_up = stack.enter_context(open(up, 'w')) if up else sys.stderr
        out_down = stack.enter_context(open(down, 'w')) if down else sys.stderr

        for line in fh:
            line_split = line.split()
            p_adj = float(line_split[9])
            m_adj = float(line_split[8])

            if abs(m_adj) >= log_fc and abs(p_adj) <= p_value:
                if m_adj > 0:
                    out_up.write(line)
                else:
                    out_down.write(line)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='FILE', nargs='?', default=[],
        help='Output file of HiCcompare (default: stdin)')
    custom.add_argument(
        '--p_value', default=0.05, type=float,
        help='P-value threshold for filtering  (default: %(default)s)')
    custom.add_argument(
        '--log_fc', default=0, type=float,
        help='Absolute log FC threshold for filtering (default: %(default)s)')
    custom.add_argument(
        '--up', default=None,
        help='Output file for UP interactions (default: stderr)')
    custom.add_argument(
        '--down', default=None,
        help='Output file for DOWN interactions (default: stderr)')
    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'

    base = argparse.ArgumentParser(add_help=False)
    base.add_argument(
        '--version', action='version', version=f'%(prog)s {__version__}')
    base.add_argument(
        '--verbose', action='store_const', const=logging.DEBUG,
        default=logging.INFO, help='verbose logging for debugging')

    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[base, custom])
    args = parser.parse_args()

    log_format = '%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
    logging.basicConfig(level=args.verbose, format=log_format)

    return args


if __name__ == '__main__':
    args = parse_arguments()
    return_code = args.function(**vars(args))
    logging.shutdown()
    sys.exit(return_code)
