#!/usr/bin/env python3

import os
import numpy as np
import argparse

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
        '--insulations', nargs = '*', default=None,
        help='Insulation score outputs of hicFindTADs (ending "tad_score.bm").')
    parser.add_argument(
        '--tads', nargs='*', default=[],
        help = 'TAD scores in ".links" format.')
    parser.add_argument(
        '--flip',  action='store_true',
        help='Plot a copy of the HiC map inverted.')
    parser.add_argument(
        '--loops', nargs='*', default=[],
        help = 'Loop output.')
    parser.add_argument(
        '--stripes', nargs=2,
        help = 'Forward and Reverse stripe score in bedGraph.')
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
        '--matrix2',
        help = 'If --flip argumnet is called, use to plot a different HiC'
        'matrix as inverted.')
    parser.add_argument(
        '--log_matrix2', default=False, action='store_true',
        help='Log transform counts of matrix 2.')
    parser.add_argument(
        '--bigWig', metavar='TITLE,FILE', default=[],
        type=commaPair, action='append',
        help='Add title and bigWig files as comma seperated pairs.'
        'Call multiple times to add more files.')
    parser.add_argument(
        '--sumLogFC', nargs=2,
        help='Pair of bigWig files for sum logFC of UP and DOWN interactions.')
    parser.add_argument(
        '--sumLogFC_title', default='sum(logFC)',
        help='Title for sumLogFC track')
    parser.add_argument(
        '--bed', metavar='TITLE,FILE', default=[],
        type=commaPair, action='append',
        help='Add title and bed files as comma seperated pairs.'
        'Call multiple times to add more files.')
    parser.add_argument(
        '--depth', type = int, default = 1000000,
        help = 'HiC matrix depth.')
    parser.add_argument(
        '--colourmap', default='Purples',
        help = 'Matplotlib colour map to use for the heatmap.')
    parser.add_argument(
        '--vMin', type=float,
        help = 'Minimum score value for HiC matrices.')
    parser.add_argument(
        '--vMax', type=float,
        help = 'Maximum score value for matrix.')

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


def make_config(insulations, matrix, log, matrix2, log_matrix2, tads, loops,
                bigWig, bed, compare, sumLogFC, sumLogFC_title, stripes,
                depth, colourmap, vMin, vMax, flip):


    print('[spacer]')
    if matrix and not_empty(matrix):
        write_matrix(matrix, cmap=colourmap, depth=depth,
            vMin=vMin, vMax=vMax, log=log)

    for i, loop in enumerate(loops):
        if not_empty(loop):
            write_loops(loop, i=i, compare=compare)

    for i, tad in enumerate(tads):
        if not_empty(tad):
            write_tads(tad, i = i)

    if flip:
        inverted_matrix = matrix2 if matrix2 else matrix
        log_transform = log_matrix2 if matrix2 else log
        if not_empty(inverted_matrix):
            write_matrix(inverted_matrix, cmap=colourmap, depth=depth,
                vMin=vMin, vMax=vMax, invert=True, log=log_transform)

    print('[spacer]')

    if sumLogFC is not None:
        for i, file in enumerate(sumLogFC):
            if not_empty(file):
                if i == 0:
                    overlay = 'no'
                    colour = '#FF000080'
                else:
                    overlay = 'share-y'
                    colour = '#0000FF80'
                write_bigwig(file=file, title=sumLogFC_title, type='bedgraph',
                             alpha=0.5, colour=colour, overlay=overlay)
        print('[spacer]')


    if insulations is not None:
        for i, insulation in enumerate(insulations):
            if not_empty(insulation):
                write_insulation(insulation = insulation, i = i)
        print('[spacer]')

    if stripes is not None:
        writeStripes(stripes)

    print('# End Sample Specific')
    for title, file in bigWig:
        if not_empty(file):
            if file.endswith('.bedgraph'):
                type='bedgraph'
            else:
                type='bigwig'
            write_bigwig(file=file, title=title, type=type)
        print('[spacer]')

    for title, file in bed:
        if not_empty(file):
            write_bed(file=file, title=title)
        print('[spacer]')

    print('[x-axis]')


def not_empty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


def write_matrix(
        matrix, cmap='Purples',
        depth=1000000, vMin=None, vMax=None,
        invert=False, log=False):

    if invert:
        depth = int(depth / 2)
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
    if invert:
        config.append('orientation = inverted')

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


def write_insulation(insulation, i,
        colours = ['#d95f0280', '#1b9e7780', '#d902d980', '#00000080']):
    colour = colours[i]
    overlay = 'share-y' if i > 0 else 'no'

    print(f'[Bedgraph matrix]',
          f'file = {insulation}',
          f'title = Insulation',
          f'height = 3',
          f'color = {colour}',
          f'file_type = bedgraph_matrix',
          f'type = lines',
          f'overlay_previous = {overlay}', sep = '\n')


def writeStripes(stripes):
    print(f'[Forward Stripes]',
          f'file = {stripes[0]}',
          r'title = Stripe score (+:blue, -:yellow)',
          f'color = #FFC20A80',
          f'height = 3',
          f'file_type = bedgraph',
          f'overlay_previous = no',
          f'[Reverse Stripes]',
          f'file = {stripes[1]}',
          f'color = #0C7BDC80',
          f'height = 3',
          f'file_type = bedgraph',
          f'overlay_previous = share-y',
          f'[spacer]', sep='\n')


def write_bigwig(file, title, alpha=1, colour='#33a02c', type='bigwig', overlay='no'):

    print(f'[{type} - {title}]',
          f'file = {file}',
          f'title = {title}',
          f'color = {colour}',
          f'alpha = {alpha}',
          f'height = 3',
          f'number_of_bins = 500',
          f'nans_to_zeros = True',
          f'summary_method = mean',
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


if __name__ == "__main__":
    main()
