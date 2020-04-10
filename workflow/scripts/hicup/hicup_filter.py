#!/usr/bin/env python3

""" Wrapper for hicup_filter. """

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
from general import get_filepath, move_file, set_zip, svg_file


def main():

    parser = pct.make_parser(verbose=True, version=__version__)

    parser.add_argument(
        'infile', metavar='SAM/BAM',
        help='Alignment file produced by hicup_mapper.')
    parser.add_argument(
        '--output', default=None, metavar='SAM/BAM',
        help='Filtered alignment output file (default: stdout)')
    parser.add_argument(
        '--outdir', default='.',
        help='Directory to write filtered ditags')
    parser.add_argument(
        '--summary', default=None,
        help='HiCUP filter summary file (default: stderr)')
    parser.add_argument(
        '--ditag_chart', type=svg_file, default=None,
        help='Ditag size distribution plot in SVG format')
    parser.add_argument(
        '--filter_chart', type=svg_file, default=None,
        help='Filter pie chart in SVG format')
    parser.add_argument(
        '--shortest', default=0, type=int,
        help='Minimum allowable insert size (bps) (default: 0)')
    parser.add_argument(
        '--longest', default=0, type=int,
        help='Maximum allowable insert size (bps) (default: infinite)')
    requiredNamed = parser.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--digest', required=True,
        help='Digest file produced by hicup_digester')
    parser.set_defaults(function=filter)

    return (pct.execute(parser))


def filter(infile, output, outdir, summary, ditag_chart, filter_chart,
        shortest, longest, digest):

    log = pct.create_logger()

    zip_out = set_zip(output, ext='.bam')

    with TemporaryDirectory() as tempdir:

        command = ['hicup_filter', '--shortest', f'{shortest}',
                   '--longest', f'{longest}', '--digest', digest,
                   '--outdir', tempdir, infile]
        if zip_out:
            command.insert(1, '--zip')

        subprocess.run(command, check=True, stdout=sys.stderr)

        # Glob for summary file and move
        summary_path = get_filepath(tempdir, 'hicup_filter_summary*')
        move_file(summary_path, summary, stderr=True)

        out_charts = [ditag_chart, filter_chart]
        patterns = ['*ditag_size_distribution.svg', '*filter_piechart.svg']
        for chart, pattern in zip(out_charts, patterns):
            if chart is not None:
                chart_path = get_filepath(tempdir, pattern)
                move_file(chart_path, chart)

        # Move ditat rejects directory
        reject_path = get_filepath(tempdir, 'hicup_filter_ditag_rejects*')
        move_file(reject_path, outdir)

        subprocess.run(['ls', '-l', tempdir])
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
