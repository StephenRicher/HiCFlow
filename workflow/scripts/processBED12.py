#!/usr/bin/env python3

""" Colour code genes in BED12 file in preparation for plotting """


import os
import sys
import argparse
import pandas as pd
from collections import defaultdict
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def processBED12(bed12: str, attributes: str, label: str, typeKey: str, geneID: str):

    colourMap = ({
        'protein_coding': '31,120,180',
        'lncRNA': '251,154,153',
        'processed_pseudogene': '51,160,44',
        'miRNA': '78,223,138',
        'unknown': '0,0,0'
    })
    colourMap = defaultdict(lambda: '166,206,227', colourMap)

    bed = readBED12(bed12)
    attrs = readAttrs(attributes, label, typeKey)

    bed12fields = bed.columns # save for writing later
    isCanonical = attrs['tag'].apply(
        lambda x: 'canonical' in str(x)).rename('canonical')

    bed = pd.merge(bed, isCanonical, left_on='name', right_index=True, how='left')
    bed = pd.merge(
        bed, attrs[[geneID, label, typeKey]],
        left_on='name', right_index=True, how='left')
    bed['length'] = bed['end'] - bed['start']
    bed['name'] = bed[label]
    bed['itemRgb'] = bed.apply(setColour, args=(typeKey, colourMap), axis=1)

    bed = (
        bed.sort_values(['canonical', 'length', 'blockSizes'], ascending=False)
        .groupby(geneID)
        .head(1)
        .sort_values(['chrom', 'start'])
    )
    bed[bed12fields].to_csv(sys.stdout, index=False, header=False, sep='\t')


def readBED12(file):
    dtype = ({
        'chrom': str, 'start': int, 'end': int, 'name': str,
        'score': float, 'strand': str, 'thickStart': int, 'thickEnd': int,
        'itemRgb': str, 'blockCount': int, 'blockSizes': str, 'blockStarts': str
    })
    return pd.read_csv(file, dtype=dtype, names=dtype.keys(), sep='\t')


def readAttrs(file, label, typeKey):
    dtype = {'name': str, 'key': str, 'value': str}
    attrs = pd.read_csv(file, dtype=dtype, names=dtype.keys(), sep='\t')
    attrs = attrs.pivot(index='name', columns='key', values='value')

    if label not in attrs:
        if 'ID' in attrs:
            attrs[label] = attrs['ID']
        else:
            attrs[label] = ''
    if typeKey not in attrs:
        attrs[typeKey] = attrs['unknown']
    if 'tag' not in attrs:
        attrs['tag'] = ''
    return attrs


def selectCanonical(x):
    if x['canonical'].sum() > 0:
        return x.loc[x['canonical']]
    else:
        return bed.head().sort_values(['length', 'blockCount'], ascending=False).head(1)


def setColour(x, typeKey, colourMap):
    return colourMap[x[typeKey]]


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'bed12', help='Input genes in BED12 format.')
    parser.add_argument(
        'attributes', help='Output of attrsOut from gff3ToGenePred.')
    parser.add_argument(
        '--typeKey', default='gene_type',
        help='Attributes key describing the gene type (default: %(default)s)')
    parser.add_argument(
        '--label', default='gene_name',
        help='Attributes key to label plot (default: %(default)s)')
    parser.add_argument(
        '--geneID', default='gene_id',
        help='Attributes key describing gene ID (default: %(default)s)')
    parser.set_defaults(function=processBED12)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
