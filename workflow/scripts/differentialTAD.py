#!/usr/bin/env python3

""" Compare TADs between 2 HiC matrices via resampling """

import sys
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
from collections import defaultdict
from hicmatrix import HiCMatrix as hm
from statsmodels.stats.multitest import fdrcorrection
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def resampleCompare(
    m1: str, m2: str, referenceTADs: str,
    minBias: float, alpha: float, out: str):

    hic1 = hm.hiCMatrix(m1)
    hic2 = hm.hiCMatrix(m2)
    chrom = hic1.getChrNames()[0]
    assert chrom  == hic2.getChrNames()[0]

    # Read TADs and store indexs
    names = {'chrom': str, 'start': int, 'end': int}
    TADs = pd.read_csv(referenceTADs,
        names=names.keys(), dtype=names, usecols=[0,1,2], sep='\t')
    TADs = TADs.loc[TADs['chrom'] == chrom]
    if TADs.empty:
        Path(outDiff).touch()
        return 0

    # Read matrices and mask lower
    adjIF1 = np.triu(hic1.matrix.toarray(), k=0)
    adjIF2 = np.triu(hic2.matrix.toarray(), k=0)
    assert adjIF1.shape == adjIF2.shape

    # True only at positions on non-zero difference
    mask = (adjIF1 != 0) & (adjIF2 !=0) & (adjIF1 != adjIF2)
    # True where IF1 greater than IF2 (in valid regions only)
    positive = mask & (adjIF1 > adjIF2)
    # Process TAD domain data
    TADtrue = {}
    for i, chrom, start, end in TADs.itertuples():
        idx = hic1.getRegionBinRange(chrom, start, end)
        if idx is None:
            continue
        n = mask[idx[0]:idx[1], idx[0]:idx[1]].sum()
        nPositive = positive[idx[0]:idx[1], idx[0]:idx[1]].sum()
        nNegative = n - nPositive
        obs = nPositive - nNegative
        bias = getBias(nPositive, nNegative)
        p = stats.binom_test(nPositive, n, p=0.5, alternative='two-sided')
        TADtrue[(chrom, start, end)] = (n, nPositive, bias, obs, p)
    # Save TAD true results
    TADtrue = pd.DataFrame(TADtrue).T.reset_index()
    TADtrue.columns = ['chrom', 'start', 'end', 'n', 'totalPositive', 'bias', 'obs', 'p']
    TADtrue['totalPositive'] = TADtrue['totalPositive'].astype(int)
    TADtrue['n'] = TADtrue['n'].astype(int)
    TADtrue['p(adj)'] = fdrcorrection(TADtrue['p'])[1]
    TADtrue['group'] = TADtrue.apply(setGroup, axis=1, args=(minBias, alpha))
    TADtrue['direction'] = TADtrue.apply(setDirection, axis=1)

    cols = ['chrom', 'start', 'end', 'group', 'direction', 'p(adj)', 'p', 'bias']
    TADtrue.loc[TADtrue['group'] == 'diffTAD', cols].to_csv(
        sys.stdout, header=False, index=False, sep='\t')

    TADtrue.to_pickle(out)


def setGroup(x, minBias, pThreshold):
    if (x['bias'] > minBias) & (x['p(adj)'] < pThreshold):
        return 'diffTAD'
    else:
        return 'TAD'


def setDirection(x):
    return -1 if x['obs'] < 0 else 1


def getBias(nPos, nNeg):
    if nPos < nNeg:
        bias = 1 - (nPos / nNeg)
    else:
        bias = 1 - (nNeg / nPos)
    return bias


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=resampleCompare)
    parser.add_argument('m1', help='HiC matrix to compare in .h5 format.')
    parser.add_argument('m2', help='HiC matrix to compare in .h5 format.')
    parser.add_argument('referenceTADs', help='Reference domains to compared.')
    parser.add_argument(
        '--minBias', type=float, default=0.25,
        help='Minimum magnitude change, between 0 and 1 (default: %(default)s)')
    parser.add_argument(
        '--alpha', type=float, default=0.01,
        help='Threshold for significance (default: %(default)s)')
    parser.add_argument(
        '--out', help='Write processed data as pickle.')

    return setDefaults(parser)

    return parent


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
