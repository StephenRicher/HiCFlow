#!/usr/bin/env python3

""" Wrapper for hicup_truncater. """

__author__ = 'Stephen Richer'
__copyright__ = 'Copyright 2020, Stephen Richer'
__email__ = 'sr467@bath.ac.uk'
__license__ = 'MIT'
__version__ = '1.0.0'

import os
import sys
import logging
import argparse
import subprocess
from tempfile import TemporaryDirectory
from general import get_filepath, move_file, set_zip, restriction_seq


def main(files, output, genome, re1, re1_name, re2, re2_name, arima, **kwargs):

    if not arima and not re1:
        logging.error('Either --arima or --re1 must be set.')
        return 1
    elif arima and re1:
        logging.error('Only one of --arima and --re1 can be set.')
        return 1

    zip_out = set_zip(output, ext='.gz')

    # Write tempdir to same location as intended output to ensure enough space
    with TemporaryDirectory(dir=os.path.dirname(output)) as tempdir:

        command = ['hicup_digester',
            '--genome', genome, '--outdir', tempdir] + files
        if zip_out:
            command.insert(1, '--zip')
        if arima:
            command.insert(1, '--arima')
        else:
            command.insert(1, '--re1')
            command.insert(2, f'{re1},{re1_name}')
        if re2:
            command.insert(1, '--re2')
            command.insert(2, f'{re2},{re2_name}')
        else:
            re2_name = 'None'

        logging.info(' '.join(command))
        subprocess.run(command, check=True, stdout=sys.stderr)

        filepath = get_filepath(tempdir, 'Digest*')
        move_file(filepath, output)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'files', nargs='+', metavar='FASTA',
        help='Input reference FASTA files')
    custom.add_argument(
        '--output',
        help='Digested reference genome (default: stdout)')
    custom.add_argument(
        '--genome', default='unspecified_genome',
        help='Name of the genome to be digested')
    custom.add_argument(
        '--re1', type=restriction_seq,
        help='Restriction cut sequence with ^ to indicate cut site. '
             'e.g. Mbol = ^GATC')
    custom.add_argument(
        '--re1_name', default='re1_unspecified',
        help='Restriction enzyme 1 name')
    custom.add_argument(
        '--re2', type=restriction_seq,
        help='Specify a restriction enzyme instead of sonication '
             'to shorten di-tags.')
    custom.add_argument(
        '--re2_name', default='re2_unspecified',
        help='Restriction enzyme 2 name')
    custom.add_argument(
        '--arima', default=False, action='store_true',
        help='Set the --re1 option to that used by the Arima protocol: '
             '^GATC,DpnII:G^ANTC,Arima')

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
