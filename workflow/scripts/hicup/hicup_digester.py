#!/usr/bin/env python3

""" Wrapper for hicup_truncater. """

__author__ = 'Stephen Richer'
__copyright__ = 'Copyright 2020, Stephen Richer'
__email__ = 'sr467@bath.ac.uk'
__license__ = 'MIT'
__version__ = '1.0.0'

import re
import sys
import glob
import shutil
import argparse
import subprocess
from os import path
import pyCommonTools as pct
from tempfile import TemporaryDirectory
from timeit import default_timer as timer
from general import get_filepath, move_file, set_zip, restriction_seq


def main():

    parser = pct.make_parser(verbose=True, version=__version__)

    parser.add_argument(
        'infiles', nargs='+', metavar='FASTA',
        help='Input reference FASTA files')
    parser.add_argument(
        '--output', default=None,
        help='Digested reference genome (default: stdout)')
    parser.add_argument(
        '--summary', default='HiCUP_digester-summary.txt',
        help='HiCUP digest summary file (default: %(default)s)')
    parser.add_argument(
        '--genome', default='unspecified_genome',
        help='Name of the genome to be digested')
    parser.add_argument(
        '--re1_name', default='re1_unspecified',
        help='Restriction enzyme 1 name')
    parser.add_argument(
        '--re2', default=None, type=restriction_seq,
        help='Specify a restriction enzyme instead of sonication '
             'to shorten di-tags.')
    parser.add_argument(
        '--re2_name', default='re2_unspecified',
        help='Restriction enzyme 2 name')
    parser.add_argument(
        '--arima', default=False, action='store_true',
        help='Set the --re1 option to that used by the Arima protocol: '
             '^GATC,DpnII:G^ANTC,Arima')
    requiredNamed = parser.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--re1', required=True, type=restriction_seq,
        help='Restriction cut sequence with ^ to indicate cut site. '
             'e.g. Mbol = ^GATC')
    parser.set_defaults(function=digest)

    return (pct.execute(parser))


def digest(infiles, output, summary, genome,
        re1, re1_name, re2, re2_name, arima):

    log = pct.create_logger()

    zip_out = set_zip(output, ext='.gz')

    with TemporaryDirectory() as tempdir:

        command = ['hicup_digester', '--re1', f'{re1},{re1_name}',
            '--genome', genome, '--outdir', tempdir] + infiles
        if zip_out:
            command.insert(1, '--zip')
        if arima:
            del a[1:3] # Remove --re1 flag and parameter
            command.insert(1, '--arima')
            log.info(f'Arima mode set: overwriting --re1 {re1},{re1_name}')
        if re2 is not None:
            command.insert(1, '--re2')
            command.insert(2, f'{re2},{re2_name}')
        else:
            re2_name = 'None'

        subprocess.run(command, check=True, stdout=sys.stderr)

        filepath = get_filepath(tempdir, 'Digest*')
        move_file(filepath, output)


if __name__ == '__main__':
    log = pct.create_logger()
    start = timer()
    RC = main()
    end = timer()
    log.info(f'Total time elapsed: {end - start} seconds.')
    sys.exit(RC)
