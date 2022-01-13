#!/usr/bin/env python3

import os
import argparse
import fileinput
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
        '--insulation',
        help='Insulation score outputs of hicFindTADs (ending "tad_score.bm").')
    parser.add_argument(
        '--tads', nargs='*', default=[],
        help = 'TAD scores in ".links" format.')
    parser.add_argument(
        '--loops', nargs='*', default=[],
        help = 'Loop output.')
    parser.add_argument(
        '--links', nargs=2,
        help = 'Plot 2 links files, up in diffTAD and down in diffTAD.')
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
        '--bigWig', metavar='TITLE,FILE,SIZE', default=[],
        type=commaPair, action='append',
        help='Add title and bigWig files as comma seperated pairs.'
        'Call multiple times to add more files.')
    parser.add_argument(
        '--rgbBed', metavar='TITLE,FILE,SIZE', default=[],
        type=commaPair, action='append',
        help='Add title and rgb BED files as comma seperated pairs.'
        'Call multiple times to add more files.')
    parser.add_argument(
        '--bed', metavar='TITLE,FILE,SIZE', default=[],
        type=commaPair, action='append',
        help='Add title and bed files as comma seperated pairs.'
        'Call multiple times to add more files.')
    parser.add_argument(
        '--genes', metavar='TITLE,FILE,SIZE', default=[],
        type=commaPair, action='append',
        help='Add title and bed files as comma seperated pairs.'
        'Call multiple times to add more files.')
    parser.add_argument(
        '--collapsedBed', metavar='TITLE,FILE', default=[],
        type=commaPair, action='append',
        help='Add unlablled non-overlapping BED intervals.')
    parser.add_argument(
        '--switchScore', help='BED file of Cscore switch probability')
    parser.add_argument('--CScore', help='BED file of Cscore')
    parser.add_argument('--changeScore', help='BED file of changeScore.')
    parser.add_argument(
        '--SNPdensity', help='BED file of SNP density per interval')
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
    parser.add_argument(
        '--tmpLinks', default='.tmp.merged.links',
        help = 'Temporary links file output to fix issues '
               'with loop plotting overlay.')

    args = parser.parse_args()
    func = args.function
    args_dict = vars(args)
    [args_dict.pop(key) for key in ['function']]

    return func(**vars(args))


def commaPair(value):
    ''' Split command seperated pair and return as tuple '''

    value = value.split(',')
    return tuple(v.strip() for v in value)


def make_config(insulation, matrix, log, tads, loops, SNPdensity,
                bigWig, bed, collapsedBed, compare, rgbBed,
                depth, colourmap, vMin, vMax, switchScore,
                genes, plain, vLines, links, CScore, tmpLinks,
                changeScore):

    if plain:
        loops = []
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
        if i > 1:
            break
        if notEmpty(tad):
            if compare and len(tads) == 1:
                colour = '#000000'
            else:
                colour = ['#FF000080', '#0000FF80'][i]
            writeTADs(tad, colour)

    if notEmpty(changeScore):
        writeColourBed(
            changeScore, title='Change Score', minV=-2, maxV=2, cmap='bwr')
        print('[spacer]')

    miniTrack = False
    for title, file, size in rgbBed:
        if notEmpty(file):
            writeRGBBed(file, title, size)
        miniTrack = True
        print('[spacer]')

    if notEmpty(SNPdensity):
        writeColourBed(
            SNPdensity, title='SNP Density', minV=0, maxV=1, cmap='binary')
        miniTrack = True
        print('[spacer]')

    if notEmpty(CScore):
        writeColourBed(CScore, title='Cscore', minV=-1, maxV=1, cmap='bwr')
        miniTrack = True
        print('[spacer]')

    if notEmpty(switchScore):
        writeColourBed(switchScore, title='Switch Score', minV=0, maxV=1, cmap='binary')
        miniTrack = True
        print('[spacer]')

    if not miniTrack:
        print('[spacer]')

    if notEmpty(insulation):
        writeInsulation(insulation)
        print('[spacer]')

    print('# End Sample Specific')

    if links is not None:
        with open(tmpLinks, 'w') as fout, fileinput.input(links) as fin:
            for line in fin:
                fout.write(line)
        writeLinks(tmpLinks, 0)
        for i, link in enumerate(links, 1):
            if notEmpty(link):
                writeLinks(link, i)
        print('[spacer]')

    for title, file, size in bigWig:
        if notEmpty(file):
            if file.endswith('.bedgraph'):
                type='bedgraph'
            else:
                type='bigwig'
            write_bigwig(file=file, title=title, type=type, size=size)
        print('[spacer]')


    for title, file in collapsedBed:
        if notEmpty(file):
            writeCollapsedBed(file, title)
        print('[spacer]')

    for title, file, size in genes:
        if notEmpty(file):
            writeGenes(file=file, title=title, size=size)
        print('[spacer]')

    for title, file, size in bed:
        if notEmpty(file):
            write_bed(file=file, title=title, size=size)
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


def writeLinks(links, i):
    # Hack fix to reverse print the loops again
    colours = {0: 'white', 1: 'Reds', 2: 'Blues'}
    alpha = 0 if i == 0 else 1
    overlay = 'no' if i == 0 else 'share-y'
    print(f'[Links]',
          f'file = {links}', sep = '\n')
    if i == 0:
        print('title = Loops')
    print(f'links_type = arcs',
          f'line_style = solid',
          f'color = {colours[i]}',
          f'height = 3',
          f'alpha = {alpha}',
          f'compact_arcs_level = 0',
          f'file_type = links',
          f'overlay_previous = {overlay}', sep = '\n')


def writeTADs(tads, colour):
    print(f'[Tads]',
          f'file = {tads}',
          f'file_type = bed',
          f'display = triangles',
          f'border_color = {colour}',
          f'color = none',
          f'overlay_previous = share-y',
          f'line_width = 1', sep = '\n')


def writeInsulation(insulation):
    print(f'[Bedgraph matrix]',
          f'file = {insulation}',
          f'title = Insulation',
          f'height = 3',
          f'file_type = bedgraph_matrix',
          f'type = lines',
          f'overlay_previous = no', sep = '\n')


def write_bigwig(
        file, title, size, alpha=1, colour='#33a02c',
        type='bigwig', overlay='no'):

    print(f'[{type} - {title}]',
          f'file = {file}',
          f'title = {title}',
          f'color = {colour}',
          f'alpha = {alpha}',
          f'height = {size}',
          f'nans_to_zeros = True',
          f'show_data_range = true',
          f'file_type = {type}',
          f'overlay_previous = {overlay}', sep = '\n')


def write_bed(file, title, size):
    print(f'[Bed - {title}]',
          f'file = {file}',
          f'title = {title}',
          f'height = {size}',
          f'file_type = bed',
          f'labels = true', sep = '\n')

def writeGenes(file, title, size):
    print(f'[rgb BED]',
          f'file = {file}',
          f'title = {title}',
          f'color = bed_rgb',
          f'border_color = none',
          f'height = {size}',
          f'file_type = bed',
          f'labels = true', sep='\n')


def writeVlines(bed):
    print(f'[vlines]',
          f'file = {bed}',
          f'type = vlines', sep='\n')


def writeRGBBed(bed, title, size):
    print(f'[rgb BED]',
          f'file = {bed}',
          f'title = {title}',
          f'labels = false',
          f'color = bed_rgb',
          f'border_color = none',
          f'line_width = 0',
          f'fontsize = 0',
          f'height = {size}',
          f'display = collapsed', sep='\n')


def writeColourBed(bed, title, minV=-1, maxV=1, cmap='bwr'):
    print(f'[{title}]',
          f'file = {bed}',
          f'title = {title}',
          f'labels = false',
          f'color = {cmap}',
          f'border_color = none',
          f'min_value = {minV}',
          f'max_value = {maxV}',
          f'line_width = 0',
          f'fontsize = 0',
          f'height = 1.5',
          f'display = collapsed', sep='\n')


def writeCollapsedBed(bed, title):
    print(f'[Collapsed Bed]',
          f'file = {bed}',
          f'title = {title}',
          f'labels = false',
          f'fontsize = 0',
          f'height = 1.5',
          f'display = collapsed', sep='\n')


if __name__ == "__main__":
    main()
