#!/usr/bin/env python3

""" Run Mann Whitney U on each bin of logFC matrix to determine
    if a direction is favoured. """

import sys
import argparse
import pandas as pd
from scipy import stats
from collections import defaultdict
from statsmodels.stats.multitest import fdrcorrection
from utilities import setDefaults, createMainParent, readHomer


__version__ = '1.0.0'


def directionPreferenceHiC(file: str, threshold: float):

    mat = readHomer(file, diagonal=True, sparse=False)
    chrom = mat.attrs['chrom']
    binSize = mat.attrs['binSize']

    # Get proportion of positive logFC after removing zero bins
    sparseMat = mat.loc[mat['score'] != 0].copy()
    diffChange = defaultdict(list)
    for start, df in sparseMat.groupby('start'):
        diffChange['start'].append(start)
        _, p = stats.wilcoxon(df['score'], alternative='two-sided')
        diffChange['p'].append(p)
        direction = 1 if df['score'].median() > 0 else -1
        diffChange['direction'].append(direction)
    diffChange = pd.DataFrame(diffChange)
    diffChange['fdr'] = fdrcorrection(diffChange['p'])[1]
    diffChange['chrom'] = chrom
    diffChange['end'] = diffChange['start'] + binSize
    diffChange['name'] = '.'
    diffChange = diffChange[diffChange['fdr'] < threshold]

    columns = ['chrom', 'start', 'end', 'name', 'direction']
    diffChange.reset_index()[columns].to_csv(
        sys.stdout, index=False, header=False, sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=directionPreferenceHiC)
    parser.add_argument('file', help='HiC matrix in homer format.')
    parser.add_argument(
        '--threshold', type=float, default=0.1,
        help='FDR threshold for determing if a fold change direction is '
             'preferred (default: %(default)s)')


    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
