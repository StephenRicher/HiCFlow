#!/usr/bin/env python3

""" Subset BED file by region """

import re
import sys
import argparse
import fileinput
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def subsetBED(infile: str, region: dict) -> None:

    with fileinput.input(infile) as fh:
        n = 0
        for line in fh:
            line = line.strip()
            chrom, start, end, name, score, strand = line.split('\t')
            if chrom == region['chr']:
                n += 1
                if int(start) >= region['start'] and n > 1:
                    print(prevLine)
                if int(end) > region['end']:
                    print(line)
                    break
                prevLine = line
            elif n > 0:
                # If loop didn't break before then print
                # last line of chromosome and break
                print(prevLine)
                break


def coordinates(value):
    ''' Validate input for genomic coordinates  '''

    pattern = re.compile('^[^:-]+:[0-9]+-[0-9]+$')
    if not pattern.match(value):
        raise argparse.ArgumentTypeError(
            'Expected format is CHR:START-END e.g chr1:1-1000. '
            'Chromosome name cannot contain ": -" .')
    coords = {}
    coords['chr'], coords['start'], coords['end'] = re.sub('[:-]', ' ', value).split()
    coords['start'] = int(coords['start'])
    coords['end'] = int(coords['end'])
    if not coords['start'] < coords['end']:
        raise argparse.ArgumentTypeError(
            f'Start coordinate {coords["start"]} not less '
            f'than end coordinate {coords["end"]}.')
    else:
        return coords


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=subsetBED)
    parser.add_argument('infile', nargs='?', default=[],
        help='Input BED file (default: stdin)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--region', metavar='CHR:START-END', required=True, type=coordinates,
        help='Genomic coordinates to operate on.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
