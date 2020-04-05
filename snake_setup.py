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
    if not samples['group'].str.match(r'[^-\/g]+').all():
        sys.exit(f'Invalid group definition in {samples_file}.')
    if not samples['rep'].str.match(r'\d+').all():
        sys.exit(f'Invalid replicate definition in {samples_file}.')
    if not samples['read'].str.match(r'R[12]').all():
        sys.exit(f'Invalid read definition in {samples_file}.')

    # Define single and sample names from definitions.
    samples['single'] = (samples[['group', 'rep', 'read']]
        .apply(lambda x: '-'.join(x), axis = 1))
    samples['sample'] = (samples[['group', 'rep']]
        .apply(lambda x: '-'.join(x), axis = 1))
    samples = samples.set_index(['group', 'sample', 'single'], drop = False)

    return samples


def get_grouping(samples):

    groups = {}
    for group in samples['group']:
        groups[group] = list(samples.loc[group]['rep'].unique())
    # Extract sample names
    sample_list = list(samples['sample'].unique())

    return sample_list, groups


def load_regions(regions_file):

    regions = pd.read_table(
        regions_file,
        names=['chr', 'start', 'end', 'region'],
        index_col='region',
        dtype={'start': int, 'end': int})
    regions['length'] = regions['end'] - regions['start']

    return regions
