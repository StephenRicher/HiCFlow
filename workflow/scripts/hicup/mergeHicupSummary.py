#!/usr/bin/env python3

""" Script to combine hicup summary files for MultiQC compatibility. """

__author__ = 'Stephen Richer'
__copyright__ = 'Copyright 2020, Stephen Richer'
__email__ = 'sr467@bath.ac.uk'
__license__ = 'MIT'
__version__ = '1.0.0'


import os
import sys
import logging
import argparse


def main(truncater, mapper, filter, deduplicator, **kwargs):

    write_header()
    trunc = process_summary(truncater)
    map = process_summary(mapper)

    filename = os.path.basename('-'.join(trunc[0][0].split('-')[:2]))

    if filter:
        filt = process_summary(filter)[0]
        percent_valid = (float(filt[2])/float(filt[1])) * 100
    else:
        filt = [filename, 0, 0,	0, 0, 0, 0,	0,	0, 0, 0, 0, 0]
        percent_valid = 0

    if deduplicator:
        dedup = process_summary(deduplicator)[0]
        percent_uniques = (float(dedup[2])/float(dedup[1])) * 100
        percent_unique_trans = (float(dedup[5])/float(dedup[1])) * 100
        percent_passed = (float(dedup[2])/float(trunc[0][1])) * 100
    else:
        dedup =  [filename, 0, 0, 0, 0,	0]
        percent_uniques = 0
        percent_unique_trans = 0
        percent_passed = 0

    percent_mapped = (float(map[0][10])/float(map[0][1])) * 100

    print(filename, trunc[0][1], trunc[1][1],
          trunc[0][4], trunc[1][4], trunc[0][2], trunc[1][2],
          trunc[0][6], trunc[1][6], map[0][2], map[1][2],
          map[0][4], map[1][4], map[0][6], map[1][6],
          map[0][8], map[1][8], map[0][10], map[1][10],
          filt[2], filt[6], filt[7], filt[8],
          filt[9], filt[10], filt[11], filt[12],
          dedup[2], dedup[3], dedup[4], dedup[5],
          f'{percent_mapped:.2f}', f'{percent_valid:.2f}',
          f'{percent_uniques:.2f}', f'{percent_unique_trans:.2f}',
          f'{percent_passed:.2f}', sep='\t')


def write_header():
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


def process_summary(summary):
    with open(summary) as f:
        next(f)  # Skip header
        sample1 = next(f).strip().split()
        try:
            sample2 = next(f).strip().split()
        except StopIteration:
            sample2 = None

    return sample1, sample2


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    requiredNamed = custom.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--truncater', required=True, help='Truncater summary file.')
    requiredNamed.add_argument(
        '--mapper', required=True, help='Mapper summary file.')
    requiredNamed.add_argument(
        '--filter', help='Filter summary file.')
    requiredNamed.add_argument(
        '--deduplicator', help='Deduplicator summary file.')
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
