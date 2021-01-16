#!/usr/bin/env python3

""" Normalise and statistically describe BED scores """

import sys
import json
import logging
import argparse
import pandas as pd
from typing import List
from utilities import setDefaults, createMainParent, StoreDict


__version__ = '1.0.0'


def describeBed(bed: str, meta: List):
    """ Describe BED scores and write to JSON """

    dtypes = {'chrom': 'str', 'start': 'int', 'end': 'int',
              'name': 'str', 'score': 'float'}
    bed = pd.read_csv(bed, sep='\t', names=dtypes.keys(), dtype=dtypes)
    bed['normScore'] = bed['score'] / (bed['end'] - bed['start'])
    bedStats = bed['normScore'].describe().to_dict()
    for key, value in meta.items():
        bedStats[key] = value
    json.dump(bedStats, sys.stdout)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'bed', metavar='BED', help='Input BED file.')
    parser.add_argument(
        '--meta', metavar='KEY=VALUE',
        nargs='*', action=StoreDict,
        help='Meta data for BED file.')
    parser.set_defaults(function=describeBed)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
