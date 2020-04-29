#!/usr/bin/env python3

import os
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
        '-i', '--insulations', nargs = '*', default = None,
        help='Insulation score outputs of hicFindTADs (ending "tad_score.bm").')
    parser.add_argument(
        '-t', '--tads', nargs = '*', default = None,
        help = 'TAD scores in ".links" format.')
    parser.add_argument(
        '--flip', default=False, action='store_true',
        help='Plot inverted HiC map in addition to other HiC map.')
    parser.add_argument(
        '-l', '--loops', nargs = '*', default = None,
        help = 'Loop output.')
    parser.add_argument(
        '--compare', default=False, action='store_true',
        help='Generate .ini file for HiC compare.')
    parser.add_argument(
        '-m', '--matrix',
        help = 'HiC matrix.')
    parser.add_argument(
        '--links', nargs=2, default=None,
        help = 'UP and DOWN links files showing differential interactions.')
    parser.add_argument(
        '-c', '--ctcfs', nargs = '+', default = None,
        help = 'CTCF position in bigWig format.')
    parser.add_argument(
        '-r', '--ctcf_orientation',
        help = 'CTCF orientations in BED format.')
    parser.add_argument(
        '-g', '--genes',
        help = 'Genes in BED format.')
    parser.add_argument(
        '-d', '--depth', type = int, default = 1000000,
        help = 'HiC matrix depth.')
    parser.add_argument(
        '--colourmap', default='Purples',
        help = 'Matplotlib coluor map to use for the heatmap.')
    parser.add_argument(
        '--vmin', type = int, default = 0,
        help = 'Minimum score value.')
    parser.add_argument(
        '--vmax', type = int, default = 2,
        help = 'Minimum score value.')

    args = parser.parse_args()
    func = args.function
    args_dict = vars(args)
    [args_dict.pop(key) for key in ['function']]

    return func(**vars(args))


def make_config(insulations, matrix, tads, loops, links, ctcfs, compare,
                ctcf_orientation, genes, depth, colourmap, vmin, vmax, flip):


    print('[spacer]')
    if matrix and not_empty(matrix):
        write_matrix(matrix, colourmap, depth, vmin, vmax)

    if loops is not None:
        for i, loop in enumerate(loops):
            if not_empty(loop):
                write_loops(loop, i = i, compare=compare)

    if tads is not None:
        for i, tad in enumerate(tads):
            if not_empty(tad):
                write_tads(tad, i = i)

    if flip and matrix and not_empty(matrix):
        write_matrix(matrix, colourmap, depth, vmin, vmax, flip=True)
    print('[spacer]')

    if insulations is not None:
        for i, insulation in enumerate(insulations):
            if not_empty(insulation):
                write_insulation(insulation = insulation, i = i)
        print('[spacer]')


    if links is not None:
        for i, link in enumerate(links):
            if i == 0:
                overlay = False
                direction = 'up'
            else:
                overlay = True
                direction = 'down'
            if not_empty(link):
                write_links(link, overlay=overlay, direction=direction)
        print('[spacer]')


    print('# End Sample Specific')
    if ctcfs is not None:
        for i, ctcf in enumerate(ctcfs):
            if not_empty(ctcf):
                write_ctcf(ctcf = ctcf, i = i)
        print('[spacer]')

    if ctcf_orientation is not None and not_empty(ctcf_orientation):
        write_ctcf_direction(ctcf_orientation = ctcf_orientation)
        print('[spacer]')

    if genes and not_empty(genes):
        write_genes(genes = genes)
        print('[spacer]')
    print('[x-axis]')

def not_empty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0

def write_matrix(
        matrix, colourmap='Purples',
        depth=1000000, vmin=0, vmax=2, flip=False):

    orientation = 'inverted' if flip else 'None'
    print(f'[Matrix]',
          f'file = {matrix}',
          f'depth = {depth}',
          f'min_value = {vmin}',
          f'max_value = {vmax}',
          f'colormap = {colourmap}',
          f'orientation = {orientation}',
          f'show_masked_bins = true',
          f'file_type = hic_matrix', sep = '\n')


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


def write_tads(tads, i,
        colours = ['#1b9e7780', '#d95f0280', '#d902d980', '#00000080'],
        styles = ['dashed', 'solid', 'dashed', 'solid']):
    colour = colours[i]
    style = styles[i]

    print(f'[Tads]',
          f'file = {tads}',
          f'links_type = triangles',
          f'line_style = {style}',
          f'color = {colour}',
          f'overlay_previous = share-y',
          f'line_width = 1.5', sep = '\n')


def write_insulation(insulation, i,
        colours = ['#1b9e7780', '#d95f0280', '#d902d980', '#00000080']):
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


def write_links(link, direction, overlay=False):

    overlay = 'share-y' if overlay else 'no'
    colour = 'Reds' if direction == 'up' else 'Blues'
    group1, group2 = get_links_groups(link)
    if direction == 'up':
        title = f'Differential interactions - Red (UP in {group2})'
    else:
        title = ''

    print(f'[{group1} vs {group2} - DI {direction}]',
          f'file = {link}',
          f'title = {title}',
          f'line_style = dashed',
          f'color = {colour}',
          f'height = 10',
          f'overlay_previous = {overlay}',
          f'file_type = links', sep = '\n')


def get_links_groups(path):
    """ Retrieve group name pairs from links path. """
    base = os.path.basename(path).split('-')
    return base[0], base[2]


def write_ctcf(ctcf, i,
    colours = ['#FF000080', '#0000FF80', '#00FF0080' ]):
    colour = colours[i]
    overlay = 'share-y' if i > 0 else 'no'

    print(f'[CTCF - Rep {i}]',
          f'file = {ctcf}',
          f'color = {colour}',
          f'min_value = 0',
          f'max_value = 3',
          f'height = 3',
          f'number_of_bins = 500',
          f'nans_to_zeros = True',
          f'summary_method = mean',
          f'show_data_range = yes',
          f'file_type = bigwig',
          f'overlay_previous = {overlay}', sep = '\n')

def write_ctcf_direction(ctcf_orientation):
    print(f'[CTCF orientation]',
          f'file = {ctcf_orientation}',
          f'title = CTCF orientation',
          f'color = Reds',
          f'fontsize = 8',
          f'type = genes',
          f'height = 3',
          f'file_type = bed',
          f'labels = true', sep = '\n')

def write_genes(genes):
    print(f'[Genes]',
          f'file = {genes}',
          f'type = genes',
          f'height = 3',
          f'file_type = bed',
          f'labels = true', sep = '\n')


if __name__ == "__main__":
    main()
