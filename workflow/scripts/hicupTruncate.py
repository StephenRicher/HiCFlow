#!/usr/bin/env python3

""" Wrapper for hicup_truncater. """

import os
import sys
import glob
import shutil
import tempfile
import argparse
import subprocess
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def hicupTruncate(files, hicup, output, nofill, re1, tmpdir, threads):
    with tempfile.TemporaryDirectory(dir=tmpdir) as dirname:
        command = ([
            'perl', hicup, '--re1', re1, '--zip', '--threads', str(threads),
            '--outdir', dirname, files[0], files[1]
        ])
        if nofill:
            command.insert(2, '--nofill')
        subprocess.run(command)

        shutil.move(glob.glob(dirname + '/*R1.trimmed.trunc.fastq.gz')[0], output[0])
        shutil.move(glob.glob(dirname + '/*R2.trimmed.trunc.fastq.gz')[0], output[1])

        summary = glob.glob(dirname + '/hicup_truncater_summary*')[0]
        writeSummary(summary)


def writeSummary(summary):
    """ Parse truncation summary compatable with multiQC. """
    printHeader()
    trunc = processSummary(summary)
    percentMapped = 0
    map = [[0] * 11] * 2
    filename = os.path.basename('-'.join(trunc[0][0].split('-')[:2]))
    filt = [filename, 0, 0,	0, 0, 0, 0,	0,	0, 0, 0, 0, 0]
    percentValid = 0
    dedup =  [filename, 0, 0, 0, 0,	0]
    percentUnique = 0
    percentUniqueTrans = 0
    percentPassed = 0
    print(filename, trunc[0][1], trunc[1][1],
          trunc[0][4], trunc[1][4], trunc[0][2], trunc[1][2],
          trunc[0][6], trunc[1][6], map[0][2], map[1][2],
          map[0][4], map[1][4], map[0][6], map[1][6],
          map[0][8], map[1][8], map[0][10], map[1][10],
          filt[2], filt[6], filt[7], filt[8],
          filt[9], filt[10], filt[11], filt[12],
          dedup[2], dedup[3], dedup[4], dedup[5],
          f'{percentMapped:.2f}', f'{percentValid:.2f}',
          f'{percentUnique:.2f}', f'{percentUniqueTrans:.2f}',
          f'{percentPassed:.2f}', sep='\t')


def printHeader():
    sys.stdout.write(
        'File\tTotal_Reads_1\tTotal_Reads_2\tNot_Truncated_Reads_1\t'
        'Not_Truncated_Reads_2\tTruncated_Read_1\tTruncated_Read_2\t'
        'Average_Length_Truncated_1\tAverage_Length_Truncated_2\t'
        'Too_Short_To_Map_Read_1\tToo_Short_To_Map_Read_2\t'
        'Unique_Alignments_Read_1\tUnique_Alignments_Read_2\t'
        'Multiple_Alignments_Read_1\tMultiple_Alignments_Read_2\t'
        'Failed_To_Align_Read_1\tFailed_To_Align_Read_2\tPaired_Read_1\t'
        'Paired_Read_2\tValid_Pairs\tInvalid_Pairs\tSame_Circularised\t'
        'Same_Dangling_Ends\tSame_Fragment_Internal\tRe_Ligation\t'
        'Contiguous_Sequence\tWrong_Size\tDeduplication_Read_Pairs_Uniques\t'
        'Deduplication_Cis_Close_Uniques\tDeduplication_Cis_Far_Uniques\t'
        'Deduplication_Trans_Uniques\tPercentage_Mapped\tPercentage_Valid\t'
        'Percentage_Uniques\tPercentage_Unique_Trans\t'
        'Percentage_Ditags_Passed_Through_HiCUP\n')


def processSummary(summary):
    with open(summary) as f:
        next(f)  # Skip header
        sample1 = next(f).strip().split()
        try:
            sample2 = next(f).strip().split()
        except StopIteration:
            sample2 = None

    return sample1, sample2


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'files', nargs=2, metavar='FASTQ', help='Input paired FASTQ files.')
    parser.add_argument(
        '--nofill', action='store_true',
        help='Hi-C protocol did NOT include a fill-in of sticky ends '
             'prior to re-ligation and therefore reads shall be '
             'truncated at the restriction site sequence.')
    parser.add_argument(
        '--threads', default=1, type=int,
        help='Number of threads to use (default: %(default)s)')
    parser.add_argument(
        '--tmpdir', help=f'Set temporary directory (default: {tempfile.gettempdir()})')
    requiredNamed = parser.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--re1', required=True,
        help='Restriction cut sequence with ^ to indicate cut site. '
             'e.g. Mbol = ^GATC')
    requiredNamed.add_argument(
        '--output', required=True, nargs=2, metavar='FASTQ',
        help='Truncated R1 and R2 FASTQ files.')
    requiredNamed.add_argument(
        '--hicup', required=True, help='Path to hicup_truncater perl script.')
    parser.set_defaults(function=hicupTruncate)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
