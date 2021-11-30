#!/usr/bin/env python3

""" Convert GFF3 to BED with gene_name """

import sys
import argparse
import fileinput
import pandas as pd
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def GFF3toBED(GFF3: str):

    with fileinput.input(GFF3) as fh:
        for line in fh:
            # Skip header
            if line.startswith('#'):
                continue
            (seqid, source, type, start, end,
             score, strand, phased, attributes) = line.strip().split()
            if type != 'gene':
                continue
            attributes = dict(attr.split('=') for attr in attributes.split(';'))
            try:
                name = attributes['gene_id'].split('.')[0]
            except KeyError:
                name = attributes['ID']
            print(seqid, start, end, name, score, strand, sep='\t')



def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument('GFF3', nargs='?', help='GFF3 input file.')
    parser.set_defaults(function=GFF3toBED)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
