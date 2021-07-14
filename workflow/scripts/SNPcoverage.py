#!/usr/bin/env python3

""" Read VCF, compute per bin density of heterozygous SNPs """


import os
import sys
import argparse
import pandas as pd
from collections import defaultdict
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def SNPcoverage(VCF: str, binSize: int):
    coverage = defaultdict(lambda: defaultdict(int))
    with open(VCF) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            chrom, pos = line.strip().split('\t')[:2]
            binStart = getBinStart(int(pos), binSize)
            coverage[chrom][binStart] += 1
    bedgraph = {k: pd.DataFrame(v, index=[0]).T for k, v in coverage.items()}
    bedgraph = (pd.concat(bedgraph, axis=0).reset_index()
        .rename({'level_0': 'chrom', 'level_1': 'start', 0: 'count'}, axis=1))
    bedgraph['end'] = bedgraph['start'] + binSize
    bedgraph[['chrom', 'start', 'end', 'count']].to_csv(
        sys.stdout, header=False, index=False, sep='\t')


def getBinStart(pos, binSize):
    """ Round down to nearest multiple of binSize """
    return pos - (pos % binSize)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument('VCF', help='VCF file of phased homozygous SNPs')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--binSize', type=int, required=True,
        help='Binsize of corresponding matrices.')
    parser.set_defaults(function=SNPcoverage)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
