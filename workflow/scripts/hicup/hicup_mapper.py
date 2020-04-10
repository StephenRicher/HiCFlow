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
import glob
import shutil
import argparse
import subprocess
from os import path
import pyCommonTools as pct
from tempfile import TemporaryDirectory
from timeit import default_timer as timer
from general import get_filepath, move_file, set_zip, svg_file


def main():

    parser = pct.make_parser(verbose=True, version=__version__)

    parser.add_argument(
        'infiles', type=str, nargs=2, metavar='FASTQ',
        help='Input pair of FASTQ files.')
    parser.add_argument(
        '--output', default=None, metavar='SAM/BAM',
        help='Paired alignment output file (default: stdout)')
    parser.add_argument(
        '--summary', default=None,
        help='HiCUP truncation summary file (default: stderr)')
    parser.add_argument(
        '--barchart1', type=svg_file, default=None,
        help='Mapper summary barchart for FASTQ R1 in SVG format')
    parser.add_argument(
        '--barchart2', type=svg_file, default=None,
        help='Mapper summary barchart for FASTQ R2 in SVG format')
    parser.add_argument(
        '--bowtie2', default='bowtie2',
        help='Specify the path to Bowtie 2')
    parser.add_argument(
        '--format', default=None,
        choices=['Sanger', 'Solexa_Illumina_1.0',
                 'Illumina_1.3', 'Illumina_1.5'],
        help='Specify FASTQ format')
    parser.add_argument(
        '-@', '--threads', default=1, type=int,
        help='Number of threads to use (default: %(default)s)')
    requiredNamed = parser.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--index', required=True,
        help='Path to the relevant reference genome Bowtie2 index')
    parser.set_defaults(function=map)

    return (pct.execute(parser))


def map(infiles, output, index, summary, barchart1,
        barchart2, bowtie2, format, threads):


    log = pct.create_logger()

    fastq_r1 = infiles[0]
    fastq_r2 = infiles[1]
    zip_out = set_zip(output, ext='.bam')

    with TemporaryDirectory() as tempdir:

        command = ['hicup_mapper', '--index', index, '--threads', str(threads),
                   '--outdir', tempdir, '--bowtie2', bowtie2,
                   fastq_r1, fastq_r2]
        if zip_out:
            command.insert(1, '--zip')
        if format is not None:
            command.insert(1, '--format')
            command.insert(2, format)

        subprocess.run(command, check=True, stdout=sys.stderr)

        # Glob for summary file and move
        summary_path = get_filepath(tempdir, 'hicup_mapper_summary*')
        move_file(summary_path, summary, stderr=True)

        # Move barchart summary figures
        for chart, fastq in zip([barchart1, barchart2], [fastq_r1, fastq_r2]):
            if chart is not None:
                chart_path = get_filepath(tempdir, f'{fastq}.mapper_barchart.svg')
                move_file(chart_path, chart)

        # Move outputs to specified directory
        filepath = get_filepath(tempdir, '*[!svg]')
        if filepath.endswith('sam.gz') and not output.endswith('sam.gz'):
            log.error('HiCUP mapper output file not BAM compressed.'
                      'Is samtools properly installed?')
            sys.exit(1)
        move_file(filepath, output)


if __name__ == '__main__':
    log = pct.create_logger()
    start = timer()
    RC = main()
    end = timer()
    log.info(f'Total time elapsed: {end - start} seconds.')
    sys.exit(RC)
