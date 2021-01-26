#!/usr/bin/env python3

""" Compute stripe interactions """

import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import zscore
from utilities import setDefaults, createMainParent, readHomer


__version__ = '1.0.0'


def computeStripes(file: str, forwardOut: str, reverseOut: str):

    mat = readHomer(file, diagonal=False, sparse=True, distanceNorm=False)
    attributes = mat.attrs
    binSize = attributes['binSize']
    mat['seperation'] = mat['start2'] - mat['start']
    mat['direction'] = np.where(mat['start2'] > mat['start'], '+', '-')

    stripeCompare = {'start': [], 'direction': [], 'size': [], 'score': []}
    for size in range(10*binSize, (200*binSize)+1, 10*binSize):
        # Get sum of interactions (per bin) up to set interaction distance
        stripeSum = (mat.loc[abs(mat['seperation']) < size, ['start', 'direction', 'score']]
                     .groupby(['start', 'direction']).agg(['median', 'count'])).reset_index()
        # Loop through each bin and direction (+ and - per bin)
        for (start, direction), dat in stripeSum.groupby(['start', 'direction']):
            score = float(dat[('score', 'median')])
            count = int(dat[('score', 'count')])
            if direction == '+':
                oppositeStart = start + size
                oppositeDirection = '-'
            else:
                oppositeStart = start - size
                oppositeDirection = '+'
            oppositeDat = (stripeSum.loc[
                (stripeSum['start'] == oppositeStart)
                & (stripeSum['direction'] == oppositeDirection)])
            if oppositeDat.empty:
                scoreDiff = np.nan
            else:
                oppositeScore = float(oppositeDat[('score', 'median')])
                oppositeCount = int(oppositeDat[('score', 'count')])
                scoreDiff = score / oppositeScore
            stripeCompare['start'].append(start)
            stripeCompare['direction'].append(direction)
            stripeCompare['size'].append(size)
            stripeCompare['score'].append(scoreDiff)
    stripeCompare = pd.DataFrame(stripeCompare)

    # Compute Z score based on ALL differences
    stripeCompare['z'] = zscore(stripeCompare['score'], nan_policy='omit')
    # For each position/direction retrieve max Z score across distances
    stripeCompare = stripeCompare.groupby(
        ['start', 'direction']).max().fillna(0).reset_index()
    stripeCompare.insert(0, 'chrom', attributes['chrom'])
    stripeCompare.insert(2, 'end', stripeCompare['start'] + binSize)

    for direction, df in stripeCompare.groupby('direction'):
        out = forwardOut if direction == '+' else reverseOut
        (df[['chrom', 'start', 'end', 'z']]
        .to_csv(out, index=False, header=False, sep='\t'))


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=computeStripes)
    parser.add_argument('file', help='HiC matrix in homer format.')
    parser.add_argument(
        'forwardOut',
        help='Output bedgraph of forward orientation stripe score.')
    parser.add_argument(
        'reverseOut',
        help='Output bedgraph of reverse orientation stripe score.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
