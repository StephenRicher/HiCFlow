#!/usr/bin/env python3

""" Apply median filter to HiC matrix in HOMER format """

import sys
import logging
import argparse
import numpy as np
import pandas as pd
from scipy.ndimage import median_filter


__version__ = '1.0.0'


def main(file, size=3, **kwargs):
    """ Read HOMER, apply median filter and write to stdout. """

    homer = readHomer(file)
    homer['data'] = median_filter(homer['data'], size=size)
    writeHomer(homer)


def readHomer(file):
    homer = {}
    with open(file) as fh:
        homer['header'] = fh.readline()
        ncols = len(homer['header'].split('\t'))
    homer['intervals'] = pd.read_csv(
        file, delimiter='\t', usecols=[0, 1], header=0)
    homer['data'] = np.loadtxt(
        file, delimiter='\t', skiprows=1, usecols=range(2, ncols))
    return homer


def writeHomer(homer):
    combined =  pd.concat([
        homer['intervals'],
        pd.DataFrame(data=homer['data'])], axis=1)
    sys.stdout.write(homer['header'])
    combined.to_csv(sys.stdout, header=False, index=False, sep='\t')


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='HOMER', help='HiC matrix in HOMER format.')
    custom.add_argument(
        '--size', type=int, default=3, help='Window size of median filter.')
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
