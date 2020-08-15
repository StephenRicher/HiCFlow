#!/usr/bin/env python3

import sys
import pandas as pd
from collections import defaultdict

class ConfigurationError(Exception):
    pass


def unset_defaults(default, outer_key=''):
    """ Return True if default or any value in default
        nested dictionary is an empty string.
    """
    out = 0
    if default == '':
        sys.stderr.write(
            f'\033[31mNo configuration provided for {outer_key} and '
             'no default available.\033[m\n')
        out += 1
    elif isinstance(default, dict):
        for key in default:
            full_key = f'{outer_key} : {key}' if outer_key else key
            out += unset_defaults(default[key], outer_key=full_key)
    return out


def set_default(config, default, key, outer_key):

    RC = 0
    if unset_defaults(default[key], outer_key=outer_key):
        RC = 1
    else:
        config[key] = default[key]
        sys.stderr.write(
            f'\033[33mNo configuration provided for {outer_key}'
            f' - setting to default: {config[key]}.\033[m\n')

    return RC


def set_config(config, default, outer_key=''):

    RC = 0
    for key in default:
        full_key = f'{outer_key} : {key}' if outer_key else key
        try:
            if isinstance(default[key], dict):
                set_config(config[key], default[key], outer_key=full_key)
            elif config[key] is None:
                RC += set_default(config, default, key, outer_key=full_key)
            else:
                sys.stderr.write(
                    f'\033[32mSetting {full_key} to: {config[key]}\033[m\n')
        except KeyError:
            RC += set_default(config, default, key, outer_key=full_key)
        except TypeError:
            RC += set_default(config, default, key, outer_key=full_key)

    if RC > 0:
        raise ConfigurationError(
            '\033[31mInvalid configuration setting.\033[m\n')

    return config



def load_samples(samples_file):

    samples = pd.read_table(
        samples_file, sep = ',', dtype = {'rep' : str})

    # Validate read file input with wildcard definitions
    if not samples['cell_type'].str.match(r'[^-\.\/]+').all():
        sys.exit(f'Invalid cell_type definition in {samples_file}.\n'
            'Cell types must not contain the following characters: - . /')
    if not samples['group'].str.match(r'[^-\.\/g]+').all():
        sys.exit(f'Invalid group definition in {samples_file}.\n'
            'Groups must not contain the following characters: - / . g')
    if not samples['rep'].str.match(r'\d+').all():
        sys.exit(f'Invalid replicate definition in {samples_file}.\n'
            'Replicates must only contain integers.')
    if not samples['read'].str.match(r'R[12]').all():
        sys.exit(f'Invalid read definition in {samples_file}.\n'
            'Reads must only be either "R1" or "R2".')

    # Define single and sample names from definitions.
    samples['single'] = (samples[['group', 'rep', 'read']]
        .apply(lambda x: '-'.join(x), axis = 1))
    samples['sample'] = (samples[['group', 'rep']]
        .apply(lambda x: '-'.join(x), axis = 1))
    # Ensure no duplicate names
    if samples['single'].duplicated().all():
        sys.exit(f'Duplicate sample name definitions in {samples_file}.\n')

    samples = samples.set_index(
        ['cell_type', 'group', 'sample', 'single'], drop = False)

    return samples


def get_grouping(samples):

    cell_types = {}
    for cell in samples['cell_type']:
        cell_types[cell] = list(samples.xs(cell, level=0)['sample'].unique())

    groups = {}
    for group in samples['group']:
        groups[group] = list(samples.xs(group, level=1)['rep'].unique())
    # Extract sample names
    sample_list = list(samples['sample'].unique())

    return sample_list, groups, cell_types


def load_regions(regions_file):

    regions = pd.read_table(
        regions_file,
        names=['chr', 'start', 'end', 'region'],
        index_col='region',
        dtype={'start': int, 'end': int})
    regions['length'] = regions['end'] - regions['start']

    # Validate read file input with wildcard definitions
    if not regions.index.to_series().str.match(r'[^-\.\/]+').all():
        sys.exit(f'Invalid region definition in {regions_file}.\n'
            'Region names must not contain the following characters: - . /')

    return regions


def load_coords(files):
    """ Read plot coordinates from plot coordinate and region BED file. """
    coords = defaultdict(list)
    for file in files:
        if not file:
            continue
        with open(file) as fh:
            for line in fh:
                chr, start, end, region = line.strip().split()
                coords[region].append(f'{chr}_{start}_{end}')
    return coords

def load_vcf_paths(phased_vcfs, samples):

    vcfs = pd.read_table(
        phased_vcfs, squeeze=True,
        names=['cell_type', 'path'],
        index_col='cell_type',sep = ',')

    if not vcfs.index.str.match(r'[^-\.\/]+').all():
        sys.exit(f'Invalid cell_type definition in {phased_vcfs}.\n'
            'Cell types must not contain the following characters: - . /')

    if len(vcfs.index.difference(samples.index.get_level_values('cell_type'))):
        sys.exit(f'Differing cell types given in {phased_vcfs} compared '
                 'to samples file.')

    return vcfs


def load_genomes(genomes_path):
    genomes = pd.read_table(
        genomes_path, squeeze=True,
        names=['cell_type', 'path'],
        index_col='cell_type', sep = ',')
    return genomes


def get_allele_groupings(samples):

    allele_groups = defaultdict(list)
    allele_samples = []
    for sample in samples:
        sample = sample.split('-')
        group = sample[0]
        rep = sample[1]
        allele_samples.extend([f'{group}_a1-{rep}', f'{group}_a2-{rep}'])
        allele_groups[f'{group}_a1'].append(rep)
        allele_groups[f'{group}_a2'].append(rep)

    return allele_groups, allele_samples
