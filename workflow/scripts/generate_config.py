#!/usr/bin/env python3

import os
import argparse
import numpy as np


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
        '--plain', default=False, action='store_true',
        help='If set ignore PCA, insulation, loops, stripes and TAD.')
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
        '--absChange',
        help='Bedgraph file for absolute sumLogFC values')
    parser.add_argument(
        '--absChange_title', default='compare score',
        help='Title for sumLogFC track')
    parser.add_argument(
        '--absChange_p', type=float, default=0.05,
        help='P-value threshold to constain axis limits of '
             'absChange plot to only show significant interactions '
             '(default: %(default)s)')
    parser.add_argument(
        '--directionScore',
        help='Bed file of bins with directional preference.')
    parser.add_argument(
        '--nBedgraphBins', type=int, default=700,
        help='Number of bins for bedgraph input')
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


def make_config(insulations, matrix, log, matrix2, log_matrix2, tads, loops,
                bigWig, bed, compare, absChange, absChange_title,
                directionScore, absChange_p, nBedgraphBins,
                stripes, depth, colourmap, vMin, vMax, flip, plain, vLines):

    if plain:
        loops = []
        stripes = None
        tads = []
        insulations = []
        vLines = None

    print('[spacer]')
    if notEmpty(matrix):
        write_matrix(matrix, cmap=colourmap, depth=depth,
            vMin=vMin, vMax=vMax, log=log)

    for i, loop in enumerate(loops):
        if notEmpty(loop):
            write_loops(loop, i=i, compare=compare)

    for i, tad in enumerate(tads):
        if notEmpty(tad):
            write_tads(tad, i = i)

    if flip:
        inverted_matrix = matrix2 if matrix2 else matrix
        log_transform = log_matrix2 if matrix2 else log
        if notEmpty(inverted_matrix):
            write_matrix(inverted_matrix, cmap=colourmap, depth=depth,
                vMin=vMin, vMax=vMax, invert=True, log=log_transform)

    print('[spacer]')

    if notEmpty(absChange):
        # Limit vmin to ensure only scores above significance are plotted
        vMin = -np.log10(absChange_p)
        write_bigwig(
            file=absChange, title=absChange_title,
            nBedgraphBins=nBedgraphBins, type='bedgraph', alpha=1,
            colour='#000000', overlay='no', vmin=vMin, vmax=20)
        print('[spacer]')

    if notEmpty(directionScore):
        writeDirectionScore(directionScore, cmap=colourmap)
    print('[spacer]')

    for i, insulation in enumerate(insulations):
        if notEmpty(insulation):
            write_insulation(insulation=insulation, compare=compare, i=i)
        if i == len(insulations) - 1:
            print('[spacer]')

    if notEmpty(stripes):
        writeStripes(stripes)

    print('# End Sample Specific')
    for title, file in bigWig:
        if notEmpty(file):
            if file.endswith('.bedgraph'):
                type='bedgraph'
            else:
                type='bigwig'
            write_bigwig(file=file, title=title,type=type)
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


def write_bigwig(
        file, title, alpha=1, colour='#33a02c',
        type='bigwig', overlay='no', nBedgraphBins=700,
        vmin=None, vmax=None):

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
    if vmin is not None:
        print(f'min_value = {vmin}')
    if vmax is not None:
        print(f'max_value = {vmax}')


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


def writeDirectionScore(bed, cmap):
    print(f'[sumLogFCBed]',
          f'file = {bed}',
          f'labels = false',
          f'color = {cmap}',
          f'line_width = 0',
          f'fontsize = 0',
          f'display = collapsed', sep='\n')

if __name__ == "__main__":
    main()
