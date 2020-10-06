#!/usr/bin/env python3

""" Create new FastQC zip directory with a new sample name for multiQC."""

import os
import sys
import logging
import zipfile
import tempfile
import argparse


__version__ = '1.0.0'


def main(zipIn, zipOut, sample, **kwargs):

    with tempfile.TemporaryDirectory() as tmpDir:

        # Create temporary file path for storing updated 'fastqc_data.txt'
        updatedData = os.path.join(tmpDir, 'fastqc_data-updated.txt')

        # Copy Zip archive to new location but extract 'fastqc_data.txt'
        with zipfile.ZipFile(zipIn, 'r') as zipInfd, zipfile.ZipFile(zipOut, 'w') as zipOutfd:
            zipOutfd.comment = zipInfd.comment
            for item in zipInfd.infolist():
                if item.filename.endswith('fastqc_data.txt'):
                    # Retrieve archive name for adding back in later.
                    arcname = item.filename
                    # Extract data to temporary directory.
                    originalData = zipInfd.extract(item.filename, tmpDir)
                else:
                    zipOutfd.writestr(item, zipInfd.read(item.filename))

        with open(originalData) as originalf, open(updatedData, 'w') as updatef:
            for line in originalf:
                if line.startswith('Filename'):
                    header, filename = line.split()
                    filename = sample
                    updatef.write(f'{header}\t{filename}\n')
                else:
                    updatef.write(line)

        # Add updated data back to the original zip path
        with zipfile.ZipFile(zipOut, mode='a', compression=zipfile.ZIP_DEFLATED) as zf:
            zf.write(updatedData, arcname)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'zipIn', help='FastQC output zip directory.')
    custom.add_argument(
        'zipOut', help='Path to updated FastQC zip directory.')
    custom.add_argument(
        'sample', help='New sample name for multiQC parsing.')
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
