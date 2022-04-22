#!/usr/bin/env python3

import os
import argparse
import fileinput
import numpy as np
from pathlib import Path

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
        '--tads', default=[],
        help = 'TAD scores in ".links" format.')
    parser.add_argument(
        '--loops', nargs='*', default=[],
        help = 'Loop output.')
    parser.add_argument(
        '--links', nargs=2, default=[],
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
        '--genes', metavar='TITLE,FILE,SIZE', default=[],
        type=commaPair, action='append',
        help='Add title and bed files as comma seperated pairs.'
        'Call multiple times to add more files.')
    parser.add_argument(
        '--CScore', nargs='*', default=[], help='BED files of Cscore')
    parser.add_argument(
        '--SNPdensity', help='BED file of SNP density per interval')
    parser.add_argument(
        '--ratioScore', help='BED file of of ratio score.')
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
    parser.add_argument(
        '--removeXAxis', default=False, action='store_true',
        help = 'Do not plot a X axis coordinate track.')
    parser.add_argument(
        '--height', type=float, default=1,
        help='Height of BED tracks.')
    parser.add_argument(
        '--geneLabelFontSize', type=int, default=10,
        help='Font size for gene labels.')

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
                bigWig, compare, rgbBed, height,
                depth, colourmap, vMin, vMax, geneLabelFontSize,
                genes, plain, vLines, links, CScore, tmpLinks,
                removeXAxis, ratioScore):

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

    if notEmpty(tads):
        writeRGBBed(tads, title='TADs', height=2, display='stacked')

    miniTrack = False
    for i, cscore in enumerate(CScore):
        if notEmpty(cscore):
            title = setCscoreTitle(i, CScore)
            writeColourBed(
                cscore, title=title, height=height,
                minV=-1, maxV=1, cmap='bwr')
            miniTrack = True

    if notEmpty(SNPdensity):
        print('[spacer]')
        writeColourBed(
            SNPdensity, title='SNP Density', height=height,
            minV=0, maxV=1, cmap='binary')
        miniTrack = True

    if notEmpty(ratioScore):
        print('[spacer]')
        writeColourBed(
            ratioScore, title='Ratio Score', height=height,
            minV=0, maxV=0.5, cmap='binary')
        miniTrack = True

    for title, file, customheight in rgbBed:
        if notEmpty(file):
            writeRGBBed(file, title, customheight)
            print('[spacer]')
        miniTrack = True

    if not miniTrack:
        print('[spacer]')

    if notEmpty(insulation):
        writeInsulation(insulation)
        print('[spacer]')

    print('# End Sample Specific')

    isLinks = False
    for link in links:
        if notEmpty(link):
            isLinks = True

    Path(tmpLinks).touch() # Force file creation
    if isLinks:
        with open(tmpLinks, 'w') as fout, fileinput.input(links) as fin:
            for line in fin:
                fout.write(line)
        writeLinks(tmpLinks, 0)
        for i, link in enumerate(links, 1):
            if notEmpty(link):
                writeLinks(link, i)
        print('[spacer]')


    for title, file, customheight in bigWig:
        if notEmpty(file):
            if file.endswith('.bedgraph'):
                type='bedgraph'
            else:
                type='bigwig'
            write_bigwig(file=file, title=title, type=type, height=customheight)
        print('[spacer]')

    for i, (title, file, customheight) in enumerate(genes):
        if notEmpty(file):
            print('[spacer]')
            writeGenes(
                file=file, title=title,
                fontSize=geneLabelFontSize, height=customheight)

    if not removeXAxis:
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
        file, title, height, alpha=1, colour='#33a02c',
        type='bigwig', overlay='no'):

    print(f'[{type} - {title}]',
          f'file = {file}',
          f'title = {title}',
          f'color = {colour}',
          f'alpha = {alpha}',
          f'height = {height}',
          f'nans_to_zeros = True',
          f'show_data_range = true',
          f'file_type = {type}',
          f'overlay_previous = {overlay}', sep = '\n')


def writeGenes(file, title, height, style='UCSC', fontSize=10):
    print(f'[rgb BED]',
          f'file = {file}',
          f'title = {title}',
          f'color = bed_rgb',
          f'color_utr = bed_rgb',
          f'color_backbone = bed_rgb',
          f'border_color = bed_rgb',
          f'style = {style}',
          f'max_labels = 50',
          f'arrow_interval = 10',
          f'fontsize = {fontSize}',
          f'height = {height}',
          f'file_type = bed',
          f'labels = true', sep='\n')


def writeVlines(bed):
    print(f'[vlines]',
          f'file = {bed}',
          f'type = vlines', sep='\n')


def writeRGBBed(bed, title, height=1.5, display='collapsed'):
    width = 0 if display == 'stacked' else 0.5
    print(f'[rgb BED]',
          f'file = {bed}',
          f'title = {title}',
          f'labels = false',
          f'color = bed_rgb',
          f'border_color = black',
          f'line_width = {width}',
          f'fontsize = 0',
          f'height = {height}',
          f'display = {display}', sep='\n')


def writeColourBed(bed, title, height=1.5, minV=-1, maxV=1, cmap='bwr', border_color='none'):
    print(f'[{title}]',
          f'file = {bed}',
          f'title = {title}',
          f'labels = false',
          f'color = {cmap}',
          f'border_color = {border_color}',
          f'min_value = {minV}',
          f'max_value = {maxV}',
          f'line_width = 0',
          f'fontsize = 0',
          f'height = {height}',
          f'display = collapsed', sep='\n')


def writeCollapsedBed(bed, title=''):
    print(f'[Collapsed Bed]',
          f'file = {bed}',
          f'title = {title}',
          f'labels = false',
          f'color = lightgray',
          f'fontsize = 0',
          f'height = 1.5',
          f'line_width = 3',
          f'display = collapsed', sep='\n')

def setCscoreTitle(i, cscore):
    if len(cscore) == 1:
        return 'Cscore'
    elif len(cscore) == 2:
        if i == 0:
            return 'Cscore (A1)'
        else:
            return 'Cscore (A2)'
    else:
        if i == 0:
            return 'Cscore (Full)'
        elif i == 1:
            return 'A1'
        else:
            return 'A2'


if __name__ == "__main__":
    main()
