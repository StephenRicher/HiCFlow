#!/usr/bin/env python3

""" Aggregate HiCRep results and plot heatmap """

import sys
import argparse
import pandas as pd
import seaborn as sns
from typing import List
import matplotlib.pyplot as plt
from utilities import setDefaults

__version__ = '1.0.0'


def main(files: List, out: str, dpi: int, vmin: float, vmax: float, fontSize: float):

    # Set global matplotlib fontisze
    plt.rcParams.update({'font.size': fontSize})

    data = pd.concat([pd.read_csv(file) for file in files])

    # Create copy and add opposing comparisons to plot full heatmap
    dataMirror = data.copy()
    dataMirror.columns = ['sample2', 'sample1', 'h', 'scc']
    data = pd.concat([data, dataMirror], sort=False)

    # Pivot data to square matrix
    data = data.pivot(index='sample1', columns='sample2', values='scc')

    fig, ax = plt.subplots()
    ax = sns.heatmap(
        data, square=True, cmap='viridis', annot=True,
        vmin=vmin, vmax=vmax, linewidth=1, linecolor='black', ax=ax)

    # Ensure masked cells are not within 'bwr' colour map.
    ax.set_facecolor('xkcd:light grey')
    ax.xaxis.set_label_text('')
    ax.yaxis.set_label_text('')
    # Make labels of y axis horizontal
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    fig.tight_layout()
    fig.savefig(out, dpi=dpi, bbox_inches='tight')


def coeff(value):
    ivalue = float(value)
    if not -1 <= ivalue <= 1:
        raise argparse.ArgumentTypeError(f'{value} must be between -1 and 1.')
    return ivalue


def parseArgs():

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'files', nargs='*', help='HiCRep output files.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--out', required=True,
        help='Output filename for HiCRep heatmap.')
    parser.add_argument(
        '--fontSize', type=float, default=12,
        help='Font size for node name on circos plot (default: %(default)s)')
    parser.add_argument(
        '--dpi', type=int, default=300,
        help='Resolution for plot (default: %(default)s)')
    parser.add_argument(
        '--vmin', type=coeff, default=None,
        help='Minimum value of colour scale (default: %(default)s)')
    parser.add_argument(
        '--vmax', type=coeff, default=None,
        help='Maximum value of colour scale (default: %(default)s)')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(main(**vars(args)))
