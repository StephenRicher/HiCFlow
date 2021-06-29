#!/usr/bin/env python3

""" Randomly swap values between a pair of HiC matrices and write to SUTM """

import sys
import gzip
import argparse
import itertools
import numpy as np
import pandas as pd
from typing import List
from collections import defaultdict
from utilities import setDefaults, createMainParent, homer2Numpy, numpy2sutm

__version__ = '1.0.0'


def shuffleHomer(
        matrices: List, outM1: str, outM2: str,
        start: int, binSize: int, seed: int):

    m1 = homer2Numpy(matrices[0])
    m2 = homer2Numpy(matrices[1])
    m1Temp = m1.copy()

    length = m1.shape[0]
    nElements = length ** 2
    np.random.seed(seed)
    randomMask = np.random.choice([True, False], size=nElements).reshape(length, length)

    m1[randomMask] = m2[randomMask]
    m2[randomMask] = m1Temp[randomMask]

    numpy2sutm(m1, start, binSize, outM1)
    numpy2sutm(m2, start, binSize, outM2)



def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument('matrices', nargs=2, help='HiC matrix in homer format.')
    parser.add_argument(
        '--seed', default=None, type=int,
        help='Seed for random swap (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--outM1', required=True,
        help='Output path of matrix 2 in SUTM format.')
    requiredNamed.add_argument(
        '--outM2', required=True,
        help='Output path of matrix 1 in SUTM format.')
    requiredNamed.add_argument(
        '--start', type=int, required=True, help='Start position of matrix.')
    requiredNamed.add_argument(
        '--binSize', type=int, required=True, help='Bin size of matrix.')
    parser.set_defaults(function=shuffleHomer)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
