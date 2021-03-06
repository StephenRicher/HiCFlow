#!/usr/bin/env python3

import os
import argparse
import numpy as np
from hicmatrix import HiCMatrix as hm

def main():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'

    formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser = argparse.ArgumentParser(
        prog='generate_config',
        description=__doc__,
        formatter_class=formatter_class,
        epilog=epilog)

    parser.set_defaults(function = make_config)

    parser.add_argument(
        '--insulations', nargs = '*', default=[],
        help='Insulation score outputs of hicFindTADs (ending "tad_score.bm").')
    parser.add_argument(
        '--tads', nargs='*', default=[],
        help = 'TAD scores in ".links" format.')
    parser.add_argument(
        '--loops', nargs='*', default=[],
        help = 'Loop output.')
    parser.add_argument(
        '--matrix',
        help = 'HiC matrix.')
    parser.add_argument(
        '--log', default=False, action='store_true',
        help='Log transform counts of matrix.')
    parser.add_argument(
        '--compare', default=False, action='store_true',
        help='Generate .ini file for HiC compare.')
    parser.add_argument(
        '--plain', default=False, action='store_true',
        help='If set ignore PCA, insulation, loops and TAD.')
    parser.add_argument(
        '--bigWig', metavar='TITLE,FILE', default=[],
        type=commaPair, action='append',
        help='Add title and bigWig files as comma seperated pairs.'
        'Call multiple times to add more files.')
    parser.add_argument(
        '--changeScore_title', default='compare score',
        help='Title for sumLogFC track')
    parser.add_argument(
        '--changeScore',
        help='Bed file of bins with directional preference.')
    parser.add_argument(
        '--permuteScore',
        help='Bed file of bins with permutation score.')
    parser.add_argument(
        '--bed', metavar='TITLE,FILE', default=[],
        type=commaPair, action='append',
        help='Add title and bed files as comma seperated pairs.'
        'Call multiple times to add more files.')
    parser.add_argument(
        '--depth', type=int, default=1000000,
        help='HiC matrix depth.')
    parser.add_argument(
        '--colourmap', default='Purples',
        help = 'Matplotlib colour map to use for the heatmap.')
    parser.add_argument(
        '--vMin', type=float,
        help = 'Minimum score value for HiC matrices.')
    parser.add_argument(
        '--vMax', type=float,
        help = 'Maximum score value for matrix.')
    parser.add_argument(
        '--vLines', help = 'BED file to plot vertical lines.')

    args = parser.parse_args()
    func = args.function
    args_dict = vars(args)
    [args_dict.pop(key) for key in ['function']]

    return func(**vars(args))


def commaPair(value):
    ''' Split command seperated pair and return as dictionary '''

    value = value.split(',')
    title = value[0].strip()
    track = value[1].strip()

    return (title, track)


def make_config(insulations, matrix, log, tads, loops,
                bigWig, bed, compare, permuteScore, changeScore_title,
                changeScore, depth, colourmap, vMin, vMax, plain, vLines):

    if plain:
        loops = []
        tads = []
        insulations = []
        vLines = None

    print('[spacer]')
    if notEmpty(matrix):
        # Set custom vMin, vMax for comparison matrices
        if compare and (vMin is None) and (vMax is None):
            threshold = computeThreshold(
                matrix, increment=0.01, propValues=0.99)
            vMin = -threshold
            vMax = threshold
        write_matrix(matrix, cmap=colourmap, depth=depth,
            vMin=vMin, vMax=vMax, log=log)

    for i, loop in enumerate(loops):
        if notEmpty(loop):
            write_loops(loop, i=i, compare=compare)

    for i, tad in enumerate(tads):
        if notEmpty(tad):
            write_tads(tad, i = i)

    print('[spacer]')

    if notEmpty(changeScore):
        writeChangeScore(changeScore, title=changeScore_title)
        print('[spacer]')

    if notEmpty(permuteScore):
        writePermuteScore(permuteScore)
        print('[spacer]')

    for i, insulation in enumerate(insulations):
        if notEmpty(insulation):
            write_insulation(insulation=insulation, compare=compare, i=i)
        if i == len(insulations) - 1:
            print('[spacer]')

    print('# End Sample Specific')
    for title, file in bigWig:
        if notEmpty(file):
            if file.endswith('.bedgraph'):
                type='bedgraph'
            else:
                type='bigwig'
            write_bigwig(file=file, title=title, type=type)
        print('[spacer]')

    for title, file in bed:
        if notEmpty(file):
            write_bed(file=file, title=title)
        print('[spacer]')

    print('[x-axis]')

    if notEmpty(vLines):
        writeVlines(vLines)



def notEmpty(path):
    if not isinstance(path, str):
        return False
    else:
        return os.path.exists(path) and os.path.getsize(path) > 0


def computeThreshold(matrix, increment=0.01, propValues=0.99):
    """ Iteratively compute an appropriate vMin, vMax threshold to
        capture a proportion of non-zero elements """
    hic = hm.hiCMatrix(matrix)
    totalNonZero = len(hic.matrix.data)
    threshold = 0
    while True:
        withinThreshold = (abs(hic.matrix.data) < threshold).sum()
        if withinThreshold / totalNonZero > propValues:
            return threshold
        threshold += increment


def write_matrix(
        matrix, cmap='Purples', depth=1000000, vMin=None,
        vMax=None, log=False, height=None):

    config = ['[Matrix]',
              f'file = {matrix}',
              f'depth = {depth}',
              f'colormap = {cmap}',
              'show_masked_bins = true',
              'file_type = hic_matrix']

    # Do not log transform compare matrices (which have negative values)
    if log:
        config.append(f'transform = log1p')
    if vMin is not None:
        config.append(f'min_value = {vMin}')
    if vMax is not None:
        config.append(f'max_value = {vMax}')

    print(*config, sep='\n')


def write_loops(loops, i, compare=False):
    colours = ['Reds', 'Blues']
    colour = colours[i] if compare else '#FF000080'

    print(f'[Loops]',
          f'file = {loops}',
          f'links_type = loops',
          f'line_style = solid',
          f'color = {colour}',
          f'overlay_previous = share-y',
          f'file_type = links',
          f'line_width = 5', sep = '\n')


def write_tads(tads, i, colours = ['#FF000080', '#0000FF80']):
    colour = colours[i]

    print(f'[Tads]',
          f'file = {tads}',
          f'file_type = domains',
          f'links_type = triangles',
          f'border_color = {colour}',
          f'color = none',
          f'overlay_previous = share-y',
          f'line_width = 1', sep = '\n')


def write_insulation(insulation, compare, i,
        colours = ['#d95f0280', '#1b9e7780', '#d902d980', '#00000080']):
    colour = colours[i]
    overlay = 'share-y' if i > 0 else 'no'
    if compare:
        title = f'Insulation difference (Z), threshold = 2'
    else:
        title = 'Insulation'
    print(f'[Bedgraph matrix]',
          f'file = {insulation}',
          f'title = {title}',
          f'height = 3',
          f'color = {colour}',
          f'file_type = bedgraph_matrix',
          f'type = lines',
          f'overlay_previous = {overlay}', sep = '\n')
    if compare:
        writeHline(2)


def write_bigwig(
        file, title, alpha=1, colour='#33a02c',
        type='bigwig', overlay='no'):

    print(f'[{type} - {title}]',
          f'file = {file}',
          f'title = {title}',
          f'color = {colour}',
          f'alpha = {alpha}',
          f'height = 3',
          f'nans_to_zeros = True',
          f'show_data_range = true',
          f'file_type = {type}',
          f'overlay_previous = {overlay}', sep = '\n')


def write_bed(file, title):
    print(f'[Bed - {title}]',
          f'file = {file}',
          f'title = {title}',
          f'type = genes',
          f'height = 3',
          f'file_type = bed',
          f'labels = true', sep = '\n')


def writeHline(y):
    print(f'[Hlines overlayed]',
          f'color = black',
          f'line_style = dotted',
          f'y_values = {y}',
          f'overlay_previous = share-y',
          f'file_type = hlines', sep='\n')


def writeVlines(bed):
    print(f'[vlines]',
          f'file = {bed}',
          f'type = vlines', sep='\n')


def writeChangeScore(bed, title):
    print(f'[sumLogFCBed]',
          f'file = {bed}',
          f'title = {title}',
          f'labels = false',
          f'color = bed_rgb',
          f'border_color = none',
          f'line_width = 0',
          f'fontsize = 0',
          f'height = 1.5',
          f'display = collapsed', sep='\n')


def writePermuteScore(bed):
    print(f'[Permute Score]',
          f'file = {bed}',
          f'title = Permute score',
          f'labels = false',
          f'color = binary',
          f'border_color = none',
          f'min_value = 0',
          f'max_value = 1',
          f'line_width = 0',
          f'fontsize = 0',
          f'height = 1.5',
          f'display = collapsed', sep='\n')


if __name__ == "__main__":
    main()
