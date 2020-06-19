#!/usr/bin/env python3

""" Wrapper for hicup_truncater. """

__author__ = 'Stephen Richer'
__copyright__ = 'Copyright 2020, Stephen Richer'
__email__ = 'sr467@bath.ac.uk'
__license__ = 'MIT'
__version__ = '1.0.0'

import os
import sys
import shutil
import logging
import argparse
import subprocess
from tempfile import TemporaryDirectory
from general import get_filepath, move_file, set_zip


def main(files, output, index, summary, bowtie2, format, threads, **kwargs):

    zip_out = set_zip(output, ext='.bam')

    # Write tempdir to same location as intended output to ensure enough space
    with TemporaryDirectory(dir=os.path.dirname(output)) as tempdir:

        command = ['hicup_mapper', '--index', index, '--threads', str(threads),
                   '--outdir', tempdir, '--bowtie2', bowtie2,
                   files[0], files[1]]

        if zip_out:
            command.insert(1, '--zip')
        if format is not None:
            command.insert(1, '--format')
            command.insert(2, format)

        logging.info(' '.join(command))
        subprocess.run(command, check=True, stdout=sys.stderr)

        summary_path = get_filepath(tempdir, 'hicup_mapper_summary*')
        move_file(summary_path, summary, stderr=True)

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
        'files', nargs=2, help='Input pair of FASTQ files.')
    custom.add_argument(
        '--output',
        help='Output alignment file in BAM format (default: stdout)')
    custom.add_argument(
        '--summary',
        help='HiCUP mapper summary file (default: stderr)')
    custom.add_argument(
        '--bowtie2', default='bowtie2',
        help='Specify the path to Bowtie 2 (default: %(default)s)')
    custom.add_argument(
        '--format', default=None,
        choices=['Sanger', 'Solexa_Illumina_1.0',
                 'Illumina_1.3', 'Illumina_1.5'],
        help='Specify FASTQ format')
    custom.add_argument(
        '--threads', default=1, type=int,
        help='Number of threads to use (default: %(default)s)')
    requiredNamed = custom.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--index', required=True,
        help='Path to the relevant reference genome Bowtie2 index')
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
