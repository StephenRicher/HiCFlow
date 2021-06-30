#!/usr/bin/env python3

""" Perform permutation testing by merging BAMs from both groups and randomly
    split (keeping the relative sample sizes approximately equal). Output is
    written to SUTM format and is repeated a determined number of times. """

import os
import sys
import pysam
import random
import argparse
import pandas as pd
from typing import List
from collections import defaultdict
from utilities import setDefaults, createMainParent, getGroup, getBin

__version__ = '1.0.0'


def homer2sutm(
        BAMs: List, binSize: int, start: int, prefix1: str, prefix2: str,
        bamLogs: List, nShadow: int, seed: int):


    assert len(set([getGroup(bam) for bam in BAMs])) == 2
    random.seed(seed)

    # Create array to store shadowed counts
    split1 = [defaultdict(lambda: defaultdict(int)) for i in range(nShadow)]
    split2 = [defaultdict(lambda: defaultdict(int)) for i in range(nShadow)]

    if bamLogs is not None:
        group, probSplit = readBAMLogs(bamLogs)
    else:
        probSplit = 0.5

    for bam in BAMs:
        group = getGroup(bam)
        samfile = pysam.AlignmentFile(bam, 'rb')
        for n, record in enumerate(samfile.fetch(until_eof=True)):
            if n % 2 == 0:
                continue
            positions = sorted([
                getBin(record.next_reference_start, binSize, start),
                getBin(record.reference_start, binSize, start)])
            for i in range(nShadow):
                if random.random() > probSplit:
                    split1[i][positions[0]][positions[1]] += 1
                else:
                    split2[i][positions[0]][positions[1]] += 1

    # Ensure the correct split is written to the corect group
    if bamLogs and (group == getGroup(prefix1)):
        split1Out = prefix2
        split2Out = prefix1
    else:
        split1Out = prefix1
        split2Out = prefix2

    for i in range(nShadow):
        shadowSplit1 = [(k, k1, v1) for k, v in split1[i].items() for k1, v1 in v.items()]
        shadowSplit2 = [(k, k1, v1) for k, v in split2[i].items() for k1, v1 in v.items()]

        pd.DataFrame(shadowSplit1).sort_values([0, 1]).to_csv(
            f'{prefix1}-{i}-sutm.txt', sep=' ', header=False, index=False)
        pd.DataFrame(shadowSplit2).sort_values([0, 1]).to_csv(
            f'{prefix2}-{i}-sutm.txt', sep=' ', header=False, index=False)


def readBAMLogs(logs):
    """ Read bam logs and cout Hi-C contacts to
        determine group split ratio. """
    count = defaultdict(int)
    for log in logs:
        log = pd.read_csv(log, sep='\t')
        group = getGroup(str(log['File']))
        count[group] += int(log['Hi-C contacts'])
    assert len(count) == 2
    # Only need to return 1 groups information
    return group, count[group] / sum(count.values())


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'BAMs', nargs='*', help='BAM format alignment files.')
    parser.add_argument(
        '--nShadow', default=1, type=int,
        help='Number of shadow matrices to generate (default: %(default)s)')
    parser.add_argument(
        '--seed', default=None, type=int,
        help='Seed for random swap (default: %(default)s)')
    parser.add_argument(
        '--bamLogs', nargs='*', default=None,
        help='HiCexplorer QC information contain read counts. Used to'
             'keep size ratio of groups equal.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--binSize', type=int, required=True,
        help='Bin size of matrix.')
    requiredNamed.add_argument(
        '--start', type=int, required=True,
        help='Start position (0-based) of matrix.')
    requiredNamed.add_argument(
        '--prefix1', required=True,
        help='Path prefix of shadow matrices (group 1) in SUTM format.')
    requiredNamed.add_argument(
        '--prefix2', required=True,
        help='Path prefix of shadow matrices (group 2) in SUTM format.')
    parser.set_defaults(function=homer2sutm)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
