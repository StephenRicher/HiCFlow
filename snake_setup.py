#!/usr/bin/env python3

import sys
import pandas as pd


class ConfigurationError(Exception):
    pass


def set_config(config, default_config):

    RC = 0

    for key in default_config:
        try:
            config[key]
            sys.stderr.write(f'\033[32mSetting {key} to: {config[key]}\n')
        except KeyError:
            if default_config[key] == '':
                sys.stderr.write(
                    f'\033[31mNo configuration provided for {key} and '
                     'no default available.\n')
                RC = 1
            else:
                config[key] = default_config[key]
                sys.stderr.write(
                    f'\033[33mNo configuration provided for {key}.\n')
                sys.stderr.write(
                    f'\033[33mSetting {key} to default: {config[key]}.\n')

    if RC == 1:
        raise ConfigurationError(
            '\033[31mInvalid configuration setting.\033[m\n')

    sys.stderr.write('\033[m')
    return config


def load_samples(samples_file):

    samples = pd.read_table(
        samples_file, sep = ',', dtype = {'rep' : str})

    # Validate read file input with wildcard definitions
    if not samples['cell_type'].str.match(r'[^-\/]+').all():
        sys.exit(f'Invalid cell_type definition in {samples_file}.\n'
            'Cell types must not contain the following characters: - /')
    if not samples['group'].str.match(r'[^-\/g]+').all():
        sys.exit(f'Invalid group definition in {samples_file}.\n'
            'Groups must not contain the following characters: - / g')
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
    if not regions.index.to_series().str.match(r'[^-\/]+').all():
        sys.exit(f'Invalid region definition in {regions_file}.\n'
            'Region names must not contain the following characters: - /')

    return regions
