#!/usr/bin/env python3

import argparse
import pandas as pd


def setDefaults(parser, verbose=True, version=None):
    """ Add version and verbose arguments to parser """

    if version:
        parser.add_argument('--version', action='version',
            version=f'%(prog)s {version}')
    if verbose:
        parser.add_argument(
            '--verbose', action='store_const', const=logging.DEBUG,
            default=logging.INFO, help='verbose logging for debugging')

        args = parser.parse_args()
        logFormat='%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
        logging.basicConfig(level=args.verbose, format=logFormat)
        del args.verbose
    else:
        args = parser.parse_args()

    return args


def readHomer(matrix, binSize):
    """ Read Homer matrix format and convert to long format """

    mat = pd.read_csv(matrix, skiprows=1, header=None, sep='\t').drop(0, axis=1)
    # Split chromosome and start into 2 column DF
    positions = mat[1].str.split('-', expand=True)
    # Add end positions
    positions[2] = positions[1].astype(int) + binSize
    # Add full region position and set as index
    positions['pos'] = mat[1]
    positions =  positions.set_index('pos')
    # Set columns
    mat.columns = pd.Series('region').append(positions[1])
    mat['end'] = positions[1].values.astype(int)

    mat = mat.melt(id_vars=['region', 'end'], var_name='start', value_name='score')
    mat['start'] = mat['start'].astype(int)

    return positions, mat
