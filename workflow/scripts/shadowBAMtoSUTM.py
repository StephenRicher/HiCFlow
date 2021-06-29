#!/usr/bin/env python3

""" Radomly split SAM files and write to SUTM """

import os
import sys
import pysam
import random
import argparse
import pandas as pd
from typing import List
from collections import defaultdict
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def homer2sutm(BAMs: List, binSize: int, out1: str, out2: str, bamLogs: List, seed: int):

    random.seed(seed)

    split1 = defaultdict(lambda: defaultdict(int))
    split2 = defaultdict(lambda: defaultdict(int))

    if bamLogs is not None:
        group, probSplit = readBAMLogs(bamLogs)
    else:
        probSplit = 0.5

    for bam in BAMs:
        samfile = pysam.AlignmentFile(bam, 'rb')
        for i, record in enumerate(samfile.fetch(until_eof=True)):
            if i % 2 == 0:
                continue
            positions = sorted([
                getBin(record.next_reference_start, binSize),
                getBin(record.reference_start, binSize)])
            if random.random() > probSplit:
                split1[positions[0]][positions[1]] += 1
            else:
                split2[positions[0]][positions[1]] += 1

    split1 = [(k, k1, v1) for k, v in split1.items() for k1, v1 in v.items()]
    split2 = [(k, k1, v1) for k, v in split2.items() for k1, v1 in v.items()]

    # Ensure the correct split is written to the corect group
    if bamLogs and (group == getGroup(out1)):
        split1Out = out2
        split2Out = out1
    else:
        split1Out = out1
        split2Out = out2

    pd.DataFrame(split1).sort_values([0, 1]).to_csv(split1Out, header=False, index=False)
    pd.DataFrame(split2).sort_values([0, 1]).to_csv(split2Out, header=False, index=False)


def readBAMLogs(logs):
    """ Read bam logs to determine group split ratio """
    count = defaultdict(int)
    for log in logs:
        log = pd.read_csv(log, sep='\t')
        group = getGroup(log['File'])
        count[group] += int(log['Hi-C contacts'])
    assert len(count) == 2
    # Only need to return 1 groups information
    return group, count[group] / sum(count.values)


def getGroup(path):
    """ Retrieve group name from file path """
    basename = os.path.basename(path)
    return basename.split('-')[0]


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
    parser.add_argument(
        '--bamLogs', nargs='*', default=None,
        help='HiCexplorer QC information contain read counts. Used to'
             'keep size ratio of groups equal.')
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
