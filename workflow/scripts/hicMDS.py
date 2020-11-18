#!/usr/bin/env python3

""" Perform multidimensional scalling of HiC subtraction matrix """


import sys
import argparse
import pandas as pd
from sklearn.manifold import MDS
from utilities import setDefaults

__version__ = '1.0.0'


def hicMDS(file: str, binSize: int):
    positions, mat = readHomer(file, binSize)
    embedding = MDS(n_components=1, dissimilarity='precomputed')
    transformed = embedding.fit_transform(mat)
    bedgraph = pd.concat([positions, pd.DataFrame(transformed)], axis=1)
    # Write to stdout
    bedgraph.to_csv(sys.stdout,  sep='\t', index=False, header=False)


def readHomer(matrix, binSize):
    """ Read Homer matrix format and return positions and numpy array """
    mat = pd.read_csv(matrix, skiprows=1, header=None, sep='\t')
    # Extract chromsome and start coordinates
    positions = mat[0].str.split('-', expand=True)
    # Add end positions
    positions[2] = positions[1].astype(int) + binSize
    # Drop coordinates
    mat = mat.drop([0, 1], axis=1).to_numpy()

    return positions, mat


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('file', help='HiC matrix in homer format.')
    parser.add_argument('--binSize', type=int, help='Bin size for matrix.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(hicMDS(**vars(args)))
