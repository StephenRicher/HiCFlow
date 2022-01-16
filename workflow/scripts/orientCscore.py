#!/usr/bin/env python3

""" Orient Cscore sign using GC content """


import sys
import logging
import argparse
import pandas as pd
from scipy.stats import pearsonr
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def orientCscore(cscore: str, fasta: str) -> None:

    start, sequence = readFasta(fasta)
    names = {'chrom': str, 'start': int, 'end': int, 'score': float}
    cscore = pd.read_csv(
        cscore, names=names.keys(), dtype=names, skiprows=1, sep='\t')
    # Filter boundaries outside sequence boundary (truncated to binsize)
    nBases = len(sequence)
    cscore = cscore.loc[(cscore['start'] < nBases) & (cscore['end'] <= nBases)]
    cscore['gc'] = cscore.apply(getGC, args=(sequence, start), axis=1)
    rho, p = pearsonr(cscore['gc'], cscore['score'])
    if rho < 0:
        cscore['score'] *= -1
    logging.warning(f'GC content correlation: p = {p}, rho = {rho}.')
    if p > 0.0:
        logging.warning(
            f'P-value of correlation is not significant p = {p}.')

    cscore[['chrom', 'start', 'end', 'score']].to_csv(
        sys.stdout, header=False, index=False, sep='\t')


def readFasta(fasta):
    sequence = ''
    with open(fasta) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                start = line.strip('>').split(':')[1].split('-')[0]
                continue # ignore first header
            elif line.startswith('>'):
                break # only read first fasta sequence
            else:
                sequence += line.upper().strip()
    return int(start), sequence


def getGC(x, sequence, start):
    seqStart = x['start'] - start
    seqEnd = x['end'] - start
    sub = sequence[seqStart: seqEnd].upper()
    gc = (sub.count('G') + sub.count('C')) / len(sub)
    return gc


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=orientCscore)
    parser.add_argument(
        'cscore', help='Output Bedgraph file of Cscore.')
    parser.add_argument(
        'fasta', help='FASTA file of relevant chromosome.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
