#!/usr/bin/env python3

""" Plot average coverage per bin """


import os
import sys
import argparse
from typing import List
import matplotlib.pyplot as plt
from collections import defaultdict
from hicmatrix import HiCMatrix as hm
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def plotCoverage(matrices: List, out: str, nBins: int, dpi: int, fontSize: float):

    # Set global matplotlib fontisze
    plt.rcParams.update({'font.size': fontSize})
    fig, ax = plt.subplots(figsize=(16, 8))
    xLim = 0

    for matrix in matrices:
        sample, binSize = splitName(matrix)
        contactsPerBin = hm.hiCMatrix(matrix).matrix.sum(axis=1).A1
        # plot the cumulative histogram
        n, bins, patches = ax.hist(
            contactsPerBin, nBins, density=True,
            cumulative=-1, histtype='step', label=binSize)
        # Get appropriate xLim
        try:
            cutOff = min(bins[:-1][n < 0.1])
        except ValueError:
            cutOff = max(bins)
        xLim = max(xLim, cutOff)

    ax.set_xlim([0, xLim])
    ax.set_title(f'{sample}', loc='left')
    ax.legend(loc='upper right')
    ax.set_xlabel('Contacts per bin')
    ax.set_ylabel('Proportion')
    fig.tight_layout()
    fig.savefig(out, dpi=dpi, bbox_inches='tight')


def splitName(file):
    """ Return sample and binsize """
    path, name = os.path.split(file)
    name = name.split('.')[0]
    if name.endswith('-SNPsplit'):
        name = name[:-9]
    elif name.endswith('-full'):
        name = name[:-5]
    sample, binSize = name.rsplit('-', 1)
    return sample, int(binSize)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'matrices', nargs='+', help='HiC matrices in h5 format.')
    parser.add_argument(
        '--nBins', type=int, default=10000,
        help='Number of bins for histogram (default: %(default)s)')
    parser.add_argument(
        '--fontSize', type=float, default=12,
        help='Font size for node name on circos plot (default: %(default)s)')
    parser.add_argument(
        '--dpi', type=int, default=300,
        help='Resolution for plot (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--out', required=True, help='Outplot plot name.')
    parser.set_defaults(function=plotCoverage)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
