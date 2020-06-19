#!/usr/bin/env python3

""" Wrapper for hicup_filter. """

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
        '--summary', default=None,
        help='HiCUP filter summary file (default: stderr)')
    parser.add_argument(
        '--uniques_chart', type=svg_file, default=None,
        help='Summary of read counts before and after deduplication')
    parser.add_argument(
        '--cistrans_chart', type=svg_file, default=None,
        help='Summary of relative read pair locations.')
    parser.set_defaults(function=deduplicate)

    return (pct.execute(parser))


def deduplicate(infile, output, summary, uniques_chart, cistrans_chart):

    log = pct.create_logger()

    zip_out = set_zip(output, ext='.bam')

    # Write tempdir to same location as intended output to ensure enough space
    with TemporaryDirectory(dir=os.path.dirname(output)) as tempdir:

        command = ['hicup_deduplicator', '--outdir', tempdir, infile]
        if zip_out:
            command.insert(1, '--zip')

        subprocess.run(command, check=True, stdout=sys.stderr)

        # Glob for summary file and move
        summary_path = get_filepath(tempdir, 'hicup_deduplicator_summary*')
        move_file(summary_path, summary, stderr=True)

        out_charts = [uniques_chart, cistrans_chart]
        patterns = ['*uniques_barchart.svg', '*cis_trans_piechart.svg']
        for chart, pattern in zip(out_charts, patterns):
            if chart is not None:
                chart_path = get_filepath(tempdir, pattern)
                move_file(chart_path, chart)

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
