#!/usr/bin/env python3

""" Invert specified region and return modified FASTA genome. """

import re
import sys
import logging
import argparse
import fileinput

__version__ = '1.0.0'


def main(file, region, **kwargs):
    header = 0
    in_chr = False
    region_sequence = []
    region_length = region['end'] - region['start'] + 1
    with fileinput.input(file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                # If previous sequence was processed as in_chr
                if in_chr:
                    line = '\n' + line
                    in_chr = False
                elif line.lstrip('>').split()[0] == region['chr']:
                    in_chr = True
                    count = 1
                print(line)
            else:
                if in_chr:
                    for base in line:
                        if count == region['start']:
                            region_sequence.append(base)
                            if len(region_sequence) == region_length:
                                region_sequence = ''.join(region_sequence)
                                inverted = reverseComplement(region_sequence)
                                for i in inverted:
                                    print(i, end='' if count % 60 else '\n')
                                    count += 1
                                in_region = False
                        else:
                            print(base, end='' if count % 60 else '\n')
                            count += 1
                else:
                    print(line)
        if in_chr:
            sys.stdout.write('\n')


def reverseComplement(sequence):
    bases = 'ACGTUWSMKRYBDHVNZ'
    complement_bases = 'TGCAAWSKMYRVHDBNZ'
    table = str.maketrans(bases, complement_bases)
    return sequence.upper().translate(table)[::-1]


def region(value):
    ''' Validate input for region argument '''

    pattern = re.compile(r'^[^\s]+:[0-9]+-[0-9]+$')
    if not pattern.match(value):
        raise argparse.ArgumentTypeError(
            f'Expected format is CHR:START-END e.g chr:1-100')
    value = value.split(':')
    region = {}
    region['chr'] = value[0]
    region['start'] = int(value[1].split('-')[0])
    region['end'] = int(value[1].split('-')[1])

    if not region['start'] < region['end']:
        raise argparse.ArgumentTypeError(
            f'Start coordinate {region["start"]} not less '
            f'than end coordinate {region["end"]}.')
    else:
        return region


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='FASTA', nargs='?', default=[],
        help='Input genome reference in FASTA format (default: stdin)')
    custom.add_argument(
        'region', metavar='CHR:START-END', type=region,
        help='Region to invert in 1-based coordinates.')
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
