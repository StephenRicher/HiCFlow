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
import pyCommonTools as pct
from tempfile import TemporaryDirectory
from timeit import default_timer as timer
from general import get_filepath, move_file, set_zip, svg_file, restriction_seq


def main():

    parser = pct.make_parser(verbose=True, version=__version__)

    parser.add_argument(
        'infiles', nargs=2, metavar='FASTQ',
        help='Input pair of FASTQ files.')
    parser.add_argument(
        '--summary', default=None,
        help='HiCUP truncation summary file (default: stderr)')
    parser.add_argument(
        '--barchart1', type=svg_file, default=None,
        help='Truncation summary barchart for FASTQ R1 in SVG format')
    parser.add_argument(
        '--barchart2', type=svg_file, default=None,
        help='Truncation summary barchart for FASTQ R2 in SVG format')
    parser.add_argument(
        '--fill', dest='fill', action='store_true',
        help='Hi-C protocol did include a fill-in of sticky ends '
             'prior to re-ligation.')
    parser.add_argument(
        '--nofill', dest='fill', action='store_false',
        help='Hi-C protocol did NOT include a fill-in of sticky ends '
             'prior to re-ligation and therefore reads shall be '
             'truncated at the restriction site sequence.')
    parser.add_argument(
        '-@', '--threads', default=1, type=int,
        help='Number of threads to use (default: %(default)s)')
    requiredNamed = parser.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--outputs', required=True, nargs=2, metavar='FASTQ',
        help='Output file names.')
    requiredNamed.add_argument(
        '--re1', required=True, type=restriction_seq,
        help='Restriction cut sequence with ^ to indicate cut site. '
             'e.g. Mbol = ^GATC')
    parser.set_defaults(function=truncate, fill=True)

    return (pct.execute(parser))


def truncate(infiles, outputs, summary, barchart1, barchart2,
        fill, re1, threads):

    fastq_r1 = infiles[0]
    fastq_r2 = infiles[1]

    if outputs[0].endswith('.gz') != outputs[1].endswith('.gz'):
        log.error(f'Output files {outfiles} have '
                  'different gzip compression extensions.')
        sys.exit(1)

    zip_out = set_zip(outputs[0], ext='.gz')

    # Write tempdir to same location as intended output to ensure enough space
    with TemporaryDirectory(dir=os.path.dirname(outputs[0])) as tempdir:

        command = ['hicup_truncater', '--re1', re1, '--threads', str(threads),
                   '--outdir', tempdir, fastq_r1, fastq_r2]
        if zip_out:
            command.insert(1, '--zip')
        if not fill:
            command.insert(1, '--nofill')

        subprocess.run(command, check=True)

        # Glob for summary file and move
        summary_path = get_filepath(tempdir, 'hicup_truncater_summary*')
        move_file(summary_path, summary, stderr=True)

        # Move barchart summary figures
        for chart, fastq in zip([barchart1, barchart2], [fastq_r1, fastq_r2]):
            if chart is not None:
                chart_path = get_filepath(
                    tempdir, f'{fastq}.truncation_barchart.svg')
                move_file(chart_path, chart)

        # Move outputs to specified directory
        for fastq, output in zip([fastq_r1, fastq_r2], outputs):
            output_base = truncater_basename(fastq)
            extension = '.trunc.fastq.gz' if zip_out else '.trunc.fastq'
            output_path = os.path.join(tempdir, output_base + extension)
            shutil.move(output_path, output)


def truncater_basename(file_path):
    """ Returns the hicup_truncater specific basename. """

    base = os.path.basename(file_path)

    if base.endswith('.gz'):
        base = re.sub('.gz$', '', base)
    elif base.endswith('.bz2'):
        base = re.sub('.bz2$', '', base)
    if base.endswith('.fastq'):
        base = re.sub('.fastq$', '', base)
    elif base.endswith('.fq'):
        base = re.sub('.fq$', '', base)

    return base


if __name__ == '__main__':
    log = pct.create_logger()
    start = timer()
    RC = main()
    end = timer()
    log.info(f'Total time elapsed: {end - start} seconds.')
    sys.exit(RC)
