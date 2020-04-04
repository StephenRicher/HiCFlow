#!/usr/bin/env python3

import os
import sys
import pyCommonTools as pct
from contextlib import ExitStack
from pathlib import Path

def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(verbose=True, version=__version__)

    parser.add_argument(
        '--up', metavar='OUT_UP', nargs='*', default=[],
        help='All UP interaction link files.')
    parser.add_argument(
        '--down', metavar='OUT_DOWN', nargs='*', default=[],
        help='All DOWN interaction link files.')
    parser.add_argument(
        '--configs', metavar='CONFIGS', nargs='*', default=[],
        help='All config files to merge.')
    parser.add_argument(
        '--prefix', metavar='OUTDIR', default='',
        help='All config files to merge.')
    parser.set_defaults(function=initiate)

    return (pct.execute(parser))

def get_links_groups(path):
    """ Retrieve group name pairs from links path. """
    base = os.path.basename(path).split('-')
    return base[0], base[2]

def get_config_group(path):
    """ Retrieve group name from config. """
    base = os.path.basename(path).split('-')
    return base[0]

def write_links(link_path, direction,
        out = sys.stdout, group1 = '', group2 = '', overlay = False):

    overlay = 'share-y' if overlay else 'no'
    colour = 'Reds' if direction == 'up' else 'Blues'
    if direction == 'up':
        title = f'Differential interactions - Red (UP in {group2})'
    else:
        title = ''

    out.write(f'[{group1} vs {group2} - DI {direction}]\n'
          f'file = {link_path}\n'
          f'title = {title}\n'
          f'line_style = dashed\n'
          f'color = {colour}\n'
          f'height = 10\n'
          f'overlay_previous = {overlay}\n'
          f'file_type = links\n')

    if direction == 'down':
        out.write('[spacer]\n')

def write_merged_config(config1, config2, links_paths, prefix):

    with ExitStack() as stack:

        group1 = get_config_group(config1)
        group2 = get_config_group(config2)
        outpath = f'{prefix}{group1}-vs-{group2}.ini'

        config1_fh = stack.enter_context(pct.open(config1))
        config2_fh = stack.enter_context(pct.open(config2))
        out_fh = stack.enter_context(pct.open(outpath, 'w'))
        up_path = links_paths[0]
        down_path = links_paths[1]

        for line in config1_fh:
            if line.startswith('# End Sample Specific'):
                break
            out_fh.write(line)
        for line in config2_fh:
            if line.startswith('# End Sample Specific'):
                break
            out_fh.write(line)

        write_links(up_path, 'up', out_fh, group1, group2)
        write_links(down_path, 'down', out_fh, group1, group2, overlay = True)

        for line in config1_fh:
            out_fh.write(line)

def initiate(up, down, configs, prefix):

    links = {}
    for link in up + down:
        group1, group2 = get_links_groups(link)
        if group1 not in links:
            links[group1] = {}
        if group2 not in links[group1]:
            links[group1][group2] = []
        links[group1][group2].append(link)

    for config1 in configs:
        for config2 in configs:
            group1 = get_config_group(config1)
            group2 = get_config_group(config2)
            # Ensure all permutations are created for Snakemake
            outpath = f'{prefix}{group1}-vs-{group2}.ini'
            Path(outpath).touch()
            try:
                links_paths = links[group1][group2]
            except KeyError:
                continue

            write_merged_config(config1, config2, links_paths, prefix)

if __name__ == "__main__":
    main()
