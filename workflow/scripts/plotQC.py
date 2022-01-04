#!/usr/bin/env python3

""" Plot ditag length and insert size quality control plots """

import sys
import argparse
import pandas as pd
import seaborn as sns
from typing import List
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def plotQC(files: List, insertOut: str, ditagOut: str, dpi: int, norm: bool):

    #sns.set(font_scale=1.5)
    sns.set_style('whitegrid')
    data = [pd.read_csv(file, sep='\t') for file in files]
    data = pd.concat(data).reset_index(drop=True)

    data['group'] = data['sample'].apply(lambda x: x.split('-')[0])
    data['rep'] = data['sample'].apply(lambda x: x.split('-')[1])
    data['insert_size'] = data['insert_size'].abs()
    data['ditag_length'] = data['ditag_length'].abs()

    data = data.loc[(data['insert_size'] > 0) & (data['ditag_length'] > 0),]
    if norm:
        grouped = data.groupby(['group', 'rep', 'orientation'])
        # Retrieve sample size of smallest group
        smallestGroup = grouped['sample'].count().min()
        # Downsample data
        data = grouped.sample(n=smallestGroup)

    insertSize = sns.displot(
        data[data['interaction_type'] == 'cis'],
        x='insert_size', hue='orientation',
        col='rep', row='group', kind='kde',
        log_scale=True, facet_kws={'sharey' : 'row'})
    insertSize.set_axis_labels('Insert Size (bp)', 'Density')
    insertSize.tight_layout()
    insertSize.savefig(insertOut, dpi=dpi, bbox_inches='tight')

    ditagDist = sns.displot(data,
        x='ditag_length', hue='sample', kind='kde', log_scale=True)
    ditagDist.set_axis_labels('Ditag Length (bp)', 'Density')
    ditagDist.tight_layout()
    ditagDist.savefig(ditagOut, dpi=dpi, bbox_inches='tight')




def parseArgs():

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=plotQC)
    parser.add_argument(
        'files', nargs='*', help='Process HiC stats files.')
    parser.add_argument(
        '--norm', action='store_true',
        help='Downsample all groups to same value (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--insertOut', required=True,
        help='Output filename for insert size distribution.')
    requiredNamed.add_argument(
        '--ditagOut', required=True,
        help='Output filename for ditag size distribution.')
    parser.add_argument(
        '--dpi', type=int, default=300,
        help='Resolution for plot (default: %(default)s)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
