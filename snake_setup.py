#!/usr/bin/env python3

import sys
import pandas as pd

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


def filterRegions(regions, binSizes, nbins=100):
    """ Remove binsizes for regions  with too few bins """

    regionBin = {}
    for region, coords in regions.iterrows():
        for bin in binSizes:
            if (coords.length / bin) > nbins:
                try:
                    regionBin[region].append(bin)
                except KeyError:
                    regionBin[region] = [bin]
    return regionBin


def load_samples(samples_file):

    samples = pd.read_table(
        samples_file, sep = ',', dtype = {'rep' : str},
        usecols=['cell_type', 'group', 'rep', 'read', 'path'])

    # Validate read file input with wildcard definitions
    if not samples['cell_type'].str.match(r'[^-\.\/]+').all():
        sys.exit(f'Invalid cell_type definition in {samples_file}.\n'
            'Cell types must not contain the following characters: - . /')
    if not samples['group'].str.match(r'[^-\.\/a]+').all():
        sys.exit(f'Invalid group definition in {samples_file}.\n'
            'Groups must not contain the following characters: - / . a')
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
        cell_types[cell].extend(list(samples.xs(cell, level=0)['group'].unique()))

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
    coords = {}
    for file in files:
        if not file:
            continue
        with open(file) as fh:
            for line in fh:
                chr, start, end, region = line.strip().split()
                if region not in coords:
                    coords[region] = []
                coords[region].append(f'{chr}_{start}_{end}')
    return coords


def get_allele_groupings(samples):

    allele_groups = {}
    allele_samples = []
    for sample in samples:
        if sample.count('-') == 1:
            continue
        sample = sample.split('-')
        group = sample[0]
        rep = sample[1]
        allele_samples.extend([f'{group}_a1-{rep}', f'{group}_a2-{rep}'])
        if f'{group}_a1' not in allele_groups:
            allele_groups[f'{group}_a1'] = []
        if f'{group}_a2' not in allele_groups:
            allele_groups[f'{group}_a2'] = []
        allele_groups[f'{group}_a1'].append(rep)
        allele_groups[f'{group}_a2'].append(rep)

    return allele_groups, allele_samples


def processRestriction(samplesFile, restrictionSeqs):
    """ Match each sample with the correct set of restriction sequences based
        on the the experiment asignment. """
    samples = pd.read_csv(samplesFile, dtype = {'rep' : str})
    samples['sample'] = (samples[['group', 'rep']]
        .apply(lambda x: '-'.join(x), axis = 1))

    restrictionSeqsAdapt = {}
    for type in ['group', 'sample']:
        for name in samples[type].unique():
            experiment = samples.loc[samples[type] == name, 'experiment'].to_list()[0]
            try:
                restrictionSeqsAdapt[name] = restrictionSeqs[experiment]
            except KeyError:
                sys.exit(f'No restriction sequences defined '
                         f'for experiment {experiment}.')
    return restrictionSeqsAdapt


def unpackRestrictionSeqs(restrictionSeqs):
    unpackedRestrictionSeqs = {}
    for REpair in restrictionSeqs.values():
        for name, seq in REpair.items():
            if name in unpackedRestrictionSeqs:
                if unpackedRestrictionSeqs[name] != seq:
                    sys.exit(f'Sequence of {name} must be'
                              'consistent across experiments.')
            else:
                unpackedRestrictionSeqs[name] = seq
    return unpackedRestrictionSeqs


def addRestrictionAllle(restrictionSeqs):
    """ Add allele specific sample to restriction seq dictionary """
    alleleREseqs = restrictionSeqs.copy()
    for sample, REs in restrictionSeqs.items():
        if '-' in sample:
            group, rep = sample.split('-')
            for allele in ['a1', 'a2']:
                alleleSample = f'{group}_{allele}-{rep}'
                alleleREseqs[alleleSample] = REs
        else:
            for allele in ['a1', 'a2']:
                alleleSample = f'{sample}_{allele}'
                alleleREseqs[alleleSample] = REs
    return alleleREseqs


def allele2sample(sample):
    """ Convert allele specific sample name to normal sample """
    group, rep = sample.split('-')
    return f'{group[:-3]}-{rep}'
