#!/usr/bin/env python3

""" Radomly split SAM files and write to SUTM """

import sys
import pysam
import random
import argparse
import pandas as pd
from typing import List
from collections import defaultdict
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def homer2sutm(BAMs: List, binSize: int, out1: str, out2: str, seed: int):

    random.seed(seed)

    split1 = defaultdict(lambda: defaultdict(int))
    split2 = defaultdict(lambda: defaultdict(int))

    for bam in BAMs:
        samfile = pysam.AlignmentFile(bam, 'rb')
        for i, record in enumerate(samfile.fetch(until_eof=True)):
            if i % 2 == 0:
                continue
            positions = sorted([
                getBin(record.next_reference_start, binSize),
                getBin(record.reference_start, binSize)])
            if random.random() > 0.5:
                split1[positions[0]][positions[1]] += 1
            else:
                split2[positions[0]][positions[1]] += 1

    split1 = [(k, k1, v1) for k, v in split1.items() for k1, v1 in v.items()]
    pd.DataFrame(split1).sort_values([0, 1]).to_csv(out1, header=False, index=False)

    split2 = [(k, k1, v1) for k, v in split2.items() for k1, v1 in v.items()]
    pd.DataFrame(split2).sort_values([0, 1]).to_csv(out2, header=False, index=False)



def getBin(pos: int, binSize: int):
    return pos - (pos % binSize)

def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'BAMs', nargs='*', help='BAM format alignment files.')
    parser.add_argument(
        '--seed', default=None, type=int,
        help='Seed for random swap (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--binSize', type=int, required=True, help='Bin size of matrix.')
    requiredNamed.add_argument(
        '--out1', required=True, help='Output of matrix 1 in SUTM format.')
    requiredNamed.add_argument(
        '--out2', required=True, help='Output of matrix 2 in SUTM format.')
    parser.set_defaults(function=homer2sutm)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
