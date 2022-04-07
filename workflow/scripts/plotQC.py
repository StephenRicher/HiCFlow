#!/usr/bin/env python3

""" Plot ditag length and insert size quality control plots """

import sys
import argparse
import pandas as pd
import seaborn as sns
from typing import List
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def plotQC(files: List, out: str, dpi: int):

    sns.set(font_scale = 1.1)
    sns.set_style('whitegrid')

    data = [pd.read_pickle(file) for file in files]
    data = pd.concat(data).reset_index(drop=True)

    if len(data['group'].unique()) == 1:
        col_wrap = 3
        row = None
    else:
        col_wrap = None
        row = 'group'
    insertSize = sns.displot(
        data[data['cis']],
        x='insertSize', hue='orientation',
        col='rep', row=row, kind='kde',
        col_wrap=col_wrap, common_norm=False,
        log_scale=True, facet_kws={'sharey' : 'row'})
    for ax in insertSize.axes.flat:
        ax.set_title('')
    insertSize.set_axis_labels('Insert Size (bp)', 'Density')
    insertSize.set_titles('Group = {row_name}, Rep = {col_name}', loc='left')
    insertSize.tight_layout()
    insertSize.savefig(out, dpi=dpi, bbox_inches='tight')


def parseArgs():

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=plotQC)
    parser.add_argument(
        'out', help='Output filename for insert size distribution.')
    parser.add_argument(
        'files', nargs='*', help='Process HiC stats files.')
    parser.add_argument(
        '--dpi', type=int, default=300,
        help='Resolution for plot (default: %(default)s)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
