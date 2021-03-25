#!/usr/bin/env python3

""" Plot viewpoint plot using output bedgraphs
    of HiCFlow rule 'plotCompareViewpoint' """

import sys
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from typing import List
from pathlib import Path
import matplotlib.pyplot as plt
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def plotCompareViewpoint(bedgraphs: List, out: str, dpi: int):

    allData = []
    names=['viewRef', 'viewStart', 'viewEnd', 'chrom', 'start', 'end', 'score']
    for i, file in enumerate(bedgraphs):
        pathInfo = processFilename(file)
        df = pd.read_csv(file, sep='\t', names=names)
        df['sample'] = pathInfo['sample']
        df['logFC'] = pathInfo['logFC']
        allData.append(df)
        if (i > 0) and (pathInfo['viewpoint'] != prevViewpoint):
            logging.error('Viewpoints do not match between files.')
            raise ValueError
        prevViewpoint = pathInfo['viewpoint']
    allData = pd.concat(allData)
    allData['mid'] = (allData['end'] + allData['start']) / 2
    allData['viewMid'] = (allData['viewEnd'] + allData['viewStart']) / 2
    allData['distance'] = allData['mid'] - allData['viewMid']

    viewpoint = reformatCoordinates(pathInfo['viewpoint'])
    title = f'{pathInfo["region"]} at {pathInfo["binSize"]} binsize'
    fig, ax = plt.subplots()
    sns.lineplot(x='distance', y='score', hue='sample', alpha=0.5,
                 data=allData[allData['logFC']==False], ax=ax)
    ax.set_ylabel('Interactions')
    ax.set_xlabel(f'Distance (bp) from viewpoint at {viewpoint}')
    ax.set_title(title, loc='left')
    fig.tight_layout()
    fig.savefig(out, dpi=dpi, bbox_inches='tight')

    #sns.lineplot(x='distance', y='score', hue='sample',
    #             data=allData[allData['logFC']==True], ax=ax2)
    #ax2.set_ylabel('LogFC interactions')





def processFilename(name):
    """ Extract important information from filepath """
    name = Path(name).name # Remove directory
    name = name.split('-')
    if name[3] == 'logFC':
        logFC = True
        sample = '-'.join(name[:3])
    else:
        logFC = False
        sample = name[0] if name[3] == 'adjIF1' else name[2]
    return {'sample': sample, 'region': name[4],
            'viewpoint': name[5], 'binSize': name[6],
            'logFC': logFC}


def reformatCoordinates(coords):
    """ Convert coordinate string to nice formatting """
    coords = list(coords)
    coords[coords.index('_')] = ':'
    coords[coords.index('_')] = '-'
    return ''.join(coords)


def parseArgs():

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=plotCompareViewpoint)
    parser.add_argument(
        'bedgraphs', nargs=3,
        help='Viewpoint bedgraphs, adjIF1, adjIF2 and logFC.')
    parser.add_argument(
        '--dpi', type=int, default=300,
        help='Resolution for plot (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--out', required=True,
        help='Output filename for HiCRep heatmap.')
    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
