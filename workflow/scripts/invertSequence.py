#!/usr/bin/env python3

""" Invert specified region and return modified FASTA genome. """

import re
import sys
import argparse
import pyCommonTools as pct
from timeit import default_timer as timer


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__, infile=True, in_type='FASTA',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'infile', metavar='FASTA',
        help='Input genome reference in FASTA format.')
    parser.add_argument(
        'region', metavar='CHR:START-END', type=region,
        help='Region to invert in 1-based coordinates.')
    parser.set_defaults(function=invert)

    return (pct.execute(parser))


def invert(infile, region):
    header = 0
    in_chr = False
    region_sequence = []
    region_length = region['end'] - region['start'] + 1
    with pct.open(infile) as f:
        for line in f:
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
                                for i in reversed(region_sequence):
                                    print(i, end = '' if count % 60 else '\n')
                                    count += 1
                                in_region = False
                        else:
                            print(base, end = '' if count % 60 else '\n')
                            count += 1
                else:
                    print(line)
        if in_chr:
            sys.stdout.write('\n')


def region(value):
    ''' Validate input for region argument '''

    pattern = re.compile('^[^\s]+:[0-9]+-[0-9]+$')
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


if __name__ == '__main__':
    log = pct.create_logger()
    start = timer()
    RC = main()
    end = timer()
    log.info(f'Total time elapsed: {end - start} seconds.')
    sys.exit(RC)
