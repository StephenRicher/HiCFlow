#!/usr/bin/env python3

""" Wrapper for hicup_filter. """

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
from general import get_filepath, move_file, set_zip


def main(file, output, outdir, summary, shortest, longest, digest, **kwargs):

    zip_out = set_zip(output, ext='.bam')

    # Write tempdir to same location as intended output to ensure enough space
    with TemporaryDirectory(dir=os.path.dirname(output)) as tempdir:

        command = ['hicup_filter', '--shortest', f'{shortest}',
                   '--longest', f'{longest}', '--digest', digest,
                   '--outdir', tempdir, file]
        if zip_out:
            command.insert(1, '--zip')

        logging.info(' '.join(command))
        subprocess.run(command, check=True, stdout=sys.stderr)

        summary_path = get_filepath(tempdir, 'hicup_filter_summary*')
        move_file(summary_path, summary, stderr=True)

        # Move ditag rejects directory
        reject_path = get_filepath(tempdir, 'hicup_filter_ditag_rejects*')
        move_file(reject_path, outdir)

        # Move outputs to specified directory
        filepath = get_filepath(tempdir, '*[!svg]')
        if filepath.endswith('sam.gz') and not output.endswith('sam.gz'):
            logging.error('HiCUP mapper output file not BAM compressed.'
                      'Is samtools properly installed?')
            return 1
        move_file(filepath, output)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='SAM/BAM',
        help='Alignment file produced by hicup_mapper.')
    custom.add_argument(
        '--output', default=None, metavar='SAM/BAM',
        help='Filtered alignment output file (default: stdout)')
    custom.add_argument(
        '--outdir', default=os.getcwd(),
        help=f'Directory to write filtered ditags (default: {os.getcwd()})')
    custom.add_argument(
        '--summary',
        help='HiCUP filter summary file (default: stderr)')
    custom.add_argument(
        '--shortest', default=0, type=int,
        help='Minimum allowable insert size (bps) (default: %(default)s)')
    custom.add_argument(
        '--longest', default=0, type=int,
        help='Maximum allowable insert size (bps) (default: None)')
    requiredNamed = custom.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--digest', required=True,
        help='Digest file produced by hicup_digester')
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
