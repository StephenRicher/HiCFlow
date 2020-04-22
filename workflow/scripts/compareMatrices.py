#!/usr/bin/env python3

import os
import argparse
import subprocess

def main():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'

    formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser = argparse.ArgumentParser(
        prog='subtractMatrices',
        description=__doc__,
        formatter_class=formatter_class,
        epilog=epilog)
    parser.set_defaults(function=subtractMatrices)
    parser.add_argument(
        'matrices', nargs='+',
        help='Name of the matrices in .h5 format to use. '
             'All pairwise comparisons will be made.')
    parser.add_argument(
        '--operation', default='log2ratio',
        choices=['diff', 'ratio', 'log2ratio'],
        help='Transformation to apply to counts (default: %(default)s)')

    args = parser.parse_args()
    func = args.function
    args_dict = vars(args)
    [args_dict.pop(key) for key in ['function']]

    return func(**vars(args))


def subtractMatrices(matrices, operation='log2ratio'):
    for matrix1 in matrices:
        dir = os.path.dirname(matrix1)
        base1 = os.path.basename(matrix1)
        group1 = base1.split('-')[0]
        suffix = '-'.join(base1.split('-')[1:])
        for matrix2 in matrices:
            base2 = os.path.basename(matrix2)
            group2 = base2.split('-')[0]
            out = f'{dir}/{group1}-vs-{group2}-{suffix}'
            command = ['hicCompareMatrices', '--matrices', matrix1, matrix2,
                       '--outFileName', out, '--operation', operation]
            subprocess.run(command)


if __name__ == "__main__":
    main()
