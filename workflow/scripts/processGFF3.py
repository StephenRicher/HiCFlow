#!/usr/bin/env python3

""" Colour code genes in GFF3 file in preparation for plotting """


import os
import sys
import gzip
import argparse
from collections import defaultdict
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def processGFF3(gff3: str, label: str, typeKey: str):

    colourMap = ({
        'protein_coding': '31,120,180',
        'lncRNA': '251,154,153',
        'processed_pseudogene': '51,160,44',
        'miRNA': '78,223,138'
    })
    colourMap = defaultdict(lambda: '166,206,227', colourMap)

    if gff3.endswith('.gz'):
        read = gzip.open
        mode = 'rt'
    else:
        read = open
        mode = 'r'

    with read(gff3, mode) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            chrom, _, id_, start, end, _, strand, _, attributes = line.split()[:9]
            if id_ != 'gene':
                continue
            attributes = dict([a.split('=') for a in attributes.split(';')])
            type_ = attributes[typeKey] if typeKey in attributes else 'unknown'
            rgb = colourMap[type_]
            if (type_ == 'protein_coding') and (label in attributes):
                name = attributes[label]
            else:
                name = ' '
            print(chrom, start, end, name, 0, strand, start, end, rgb, sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument('gff3', help='Input genes in GFF3 format.')
    parser.add_argument(
        '--typeKey', default='gene_type',
        help='Attributes key describing the gene type (default: %(default)s)')
    parser.add_argument(
        '--label', default='',
        help='Attributes key to label plot (default: %(default)s)')
    parser.set_defaults(function=processGFF3)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
