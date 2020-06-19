#!/usr/bin/env python3

""" Wrapper for hicup_truncater. """

__author__ = 'Stephen Richer'
__copyright__ = 'Copyright 2020, Stephen Richer'
__email__ = 'sr467@bath.ac.uk'
__license__ = 'MIT'
__version__ = '1.0.0'

import os
import re
import sys
import shutil
import logging
import argparse
import subprocess
from tempfile import TemporaryDirectory
from general import get_filepath, move_file, set_zip, restriction_seq


def main(files, output, summary, nofill, re1, threads, **kwargs):

    if output[0].endswith('.gz') != output[1].endswith('.gz'):
        logging.error(f'Output files {output} have '
            'different gzip compression extensions.')
        return 1

    zip_out = set_zip(output[0], ext='.gz')

    # Write tempdir to same location as intended output to ensure enough space
    with TemporaryDirectory(dir=os.path.dirname(output[0])) as tempdir:

        command = ['hicup_truncater', '--re1', re1, '--threads', str(threads),
                   '--outdir', tempdir, files[0], files[1]]

        if zip_out:
            command.insert(1, '--zip')
        if nofill:
            command.insert(1, '--nofill')

        logging.info(' '.join(command))
        subprocess.run(command, check=True)

        summary_path = get_filepath(tempdir, 'hicup_truncater_summary*')
        move_file(summary_path, summary, stderr=True)

        # Move outputs to specified directory
        for fastq_in, fastq_out in zip(files, output):
            output_base = truncater_basename(fastq_in)
            extension = '.trunc.fastq.gz' if zip_out else '.trunc.fastq'
            output_path = os.path.join(tempdir, output_base + extension)
            shutil.move(output_path, fastq_out)


def truncater_basename(path):
    """ Returns the hicup_truncater specific basename. """

    base = os.path.basename(path)

    if base.endswith('.gz'):
        base = re.sub('.gz$', '', base)
    elif base.endswith('.bz2'):
        base = re.sub('.bz2$', '', base)
    if base.endswith('.fastq'):
        base = re.sub('.fastq$', '', base)
    elif base.endswith('.fq'):
        base = re.sub('.fq$', '', base)

    return base


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'files', nargs=2, metavar='FASTQ',
        help='Input paired FASTQ files.')
    custom.add_argument(
        '--summary',
        help='HiCUP truncation summary file (default: stderr)')
    custom.add_argument(
        '--nofill', action='store_true',
        help='Hi-C protocol did NOT include a fill-in of sticky ends '
             'prior to re-ligation and therefore reads shall be '
             'truncated at the restriction site sequence.')
    custom.add_argument(
        '--threads', default=1, type=int,
        help='Number of threads to use (default: %(default)s)')
    requiredNamed = custom.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--re1', required=True, type=restriction_seq,
        help='Restriction cut sequence with ^ to indicate cut site. '
             'e.g. Mbol = ^GATC')
    requiredNamed.add_argument(
        '--output', required=True, nargs=2, metavar='FASTQ',
        help='Truncated R1 and R2 FASTQ files.')
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
