#!/usr/bin/env python3

import sys
import itertools
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
    binRegion = invertRegionBin(regionBin)
    return regionBin, binRegion


def invertRegionBin(regionBin):
    """ Generate dict of bin keys to regions """

    binRegion = {}
    for region, bins in regionBin.items():
        for bin in bins:
            try:
                binRegion[str(bin)].append(region)
            except KeyError:
                binRegion[str(bin)] = [region]
    return binRegion

class HiCSamples:

    def __init__(self, samplesFile: str, restrictionSeqs: dict, alleleSpecific: bool, sep='\s+'):
        self.alleleSpecific = alleleSpecific
        self.table = self.readSamples(samplesFile, sep=sep)
        self.experimentRestriction = restrictionSeqs

    def readSamples(self, samplesFile, sep='\s+'):
        table = pd.read_table(
            samplesFile, sep=sep, dtype={'rep': str},
            usecols=['experiment', 'cell_type', 'group', 'rep', 'R1', 'R2'])

        # Validate read file input with wildcard definitions
        if not table['cell_type'].str.match(r'[^-\.\/]+').all():
            raise ValueError(
                'Cell types must not contain the following characters: - . /')
        if not table['group'].str.match(r'[^-\.\/]+').all():
            raise ValueError(
                'Groups must not contain the following characters: - . /')
        if (table.groupby('group')['cell_type'].nunique() > 1).any():
            raise ValueError(
                'An experimental group must contain only 1 cell type.')

        table['sample'] = (table[['group', 'rep']].apply(lambda x: '-'.join(x), axis=1))

        # Ensure no duplicate names
        if table['sample'].duplicated().any():
            raise ValueError(f'Duplicate sample name definitions in {samplesFile}.')

        return table


    def cellTypes(self):
        """ Return dict mapping cellType to sample """
        return self.table.groupby('cell_type')['sample'].apply(list).to_dict()


    def originalGroups(self):
        """ Return unmodified group-rep dictionary """
        return self.table.groupby('group', sort=False)['rep'].apply(list).to_dict()


    def groups(self):
        """ Return group-rep dictionary """
        if self.alleleSpecific:
            groups = {}
            for group, reps in self.originalGroups().items():
                for alleleGroup in self.group2Allele(group):
                    groups[alleleGroup] = reps
        else:
            groups = self.originalGroups()
        return groups

    def groupName(self, name):
        """ Return true is name is a group name """
        return name.count('-') == 0


    def originalSamples(self):
        """ Return unmodified sample list """
        return self.table['sample']


    def allOriginal(self):
        """ Return all unmodified sample and group names """
        return list(self.originalSamples()) + list(self.originalGroups().keys())

    def samples(self, all=False):
        """ Return sample list """
        if self.alleleSpecific:
            samples = []
            for sample in self.originalSamples():
                samples.extend(list(self.sample2Allele(sample)))
                if all:
                    samples.append(sample)
        else:
            samples = self.originalSamples().to_list()
        return samples


    def all(self):
        """ Return groups and samples as combined list """
        return self.samples() + list(self.groups())


    def groupCompares(self):
        """ Return list of pairwise group comparison """
        pairs = itertools.combinations(list(self.groups()), 2)
        return [f'{i[0]}-vs-{i[1]}' for i in pairs]

    def sampleCompares(self):
        """ Return list of pairwise sample comparison """
        pairs = itertools.combinations(self.samples(), 2)
        return [f'{i[0]}-vs-{i[1]}' for i in pairs]


    def sample2Cell(self):
        """ Return dict mapping sample name to cell type """
        sampleCell = {}
        for originalName in self.allOriginal():
            if self.groupName(originalName):
                column = 'group'
                nameConverter = self.group2Allele
            else:
                column = 'sample'
                nameConverter = self.sample2Allele
            cellType = (self.table
                .loc[self.table[column] == originalName, 'cell_type'].to_list()[0])
            nameA1, name2A2 = nameConverter(originalName)
            # Add sample and allele specific sample names to dictionary
            for name in [originalName, nameA1, name2A2]:
                sampleCell[name] = cellType
        return sampleCell


    def sample2Allele(self, sample):
        """ Convert sample name to allelic sample names """
        group, rep = sample.split('-')
        return f'{group}_a1-{rep}', f'{group}_a2-{rep}'


    def group2Allele(self, group):
        """ Convert group name """
        return f'{group}_a1', f'{group}_a2'


    def path(self, sample, read):
        """ Return path of sample-read """
        paths = self.table.loc[self.table['sample'] == sample, list(read)]
        return paths.values.tolist()[0]


    def restrictionSeqs(self, removeCut=False, dangling=False):
        """ Return dict mapping samples to restriction seqs """

        if removeCut and dangling:
            raise ValueError('Cannot set removeCut and dangling to True.')
        reSeqs = {}
        for sample in self.originalSamples():
            # Get experiment associated with sample
            experiment = (self.table
                .loc[self.table['sample'] == sample, 'experiment'].to_list()[0])
            # Get restriction seqs associated with experiment
            seqs = self.experimentRestriction[experiment]
            if dangling:
                seqs = self.restriction2dangling(seqs)
            elif removeCut:
                seqs = self.removeCutSite(seqs)
            sampleA1, sampleA2 = self.sample2Allele(sample)
            # Add sample and allele specific sample names to dictionary
            for name in [sample, sampleA1, sampleA2]:
                reSeqs[name] = seqs
        return reSeqs

    def restriction2dangling(self, seqs):
        """ Convert restriction sequence to dangling sequence """
        dangling = {}
        for name, sequence in seqs.items():
            cutIndex = sequence.index('^')
            sequence = sequence.replace('^', '')
            sequence = sequence[cutIndex:len(sequence) - cutIndex]
            dangling[name] = sequence
        return dangling

    def removeCutSite(self, seqs):
        """ Remove caret symbol from restriction sequence """
        noCut = {}
        for name, sequence in seqs.items():
            noCut[name] = sequence.replace('^', '')
        return noCut

    def restrictionNames(self, removeCut=False):
        """ Return dict mapping restriction name to sequence """
        rSeqs = {}
        for sequences in self.experimentRestriction.values():
            if removeCut:
                sequences = self.removeCutSite(sequences)
            for name, sequence in sequences.items():
                if name in rSeqs and rSeqs[name] != sequence:
                    raise ValueError(
                        f'Sequence of {name} not consistent across experiments.')
                else:
                    rSeqs[name] = sequence
        return rSeqs


def adjustCoordinates(start, end, nbases):
    """ Adjust coordinates to a multiple of nbases """
    # Round down start to closest multiple of nbases
    start = start - (start % nbases)
    # Get amount contract end position
    adjustContract = (end - start ) % nbases
    end = end - adjustContract

    return start, end


def load_regions(regions_file, adjust=1):

    regions = pd.read_table(
        regions_file,
        names=['chr', 'start', 'end', 'region'],
        index_col='region',
        dtype={'start': int, 'end': int})
    if adjust is not None:
        regions['start'], regions['end'] = adjustCoordinates(
            regions['start'], regions['end'], adjust)
    regions['length'] = regions['end'] - regions['start']

    # Validate read file input with wildcard definitions
    if not regions.index.to_series().str.match(r'[^-\.\/]+').all():
        sys.exit(f'Invalid region definition in {regions_file}.\n'
            'Region names must not contain the following characters: - . /')

    if anyOverlapping(regions):
        sys.exit(f'Overlapping intervals detected in {regions_file}.\n')

    return regions


def rangeIntersect(r1, r2):
    """ Find intersect between ranges """
    return range(max(r1.start,r2.start), min(r1.stop,r2.stop)) or None


def anyOverlapping(allRegions):
    """ Check if any within chromosome overlaps """
    for chrom, regions in allRegions.groupby('chr'):
        coordRanges = []
        for name, row in regions.iterrows():
            coordRanges.append(range(row['start'], row['end']))
        for a, b in itertools.combinations(coordRanges, 2):
            if rangeIntersect(a, b) is not None:
                return True
    return False


def reformatRegions(regions):
    """ Reformat the load_regions df to dictionary with merged
        coordinates. Serves as valid input for plot coordinates
        and plot viewpoint. """

    coords = regions[['chr', 'start', 'end']].astype(str).agg('_'.join, axis=1).to_dict()
    coords = {name:[pos] for (name, pos) in coords.items()}
    return coords


def outOfRange(chr, start, end, region):
    """ Ensure coordinates are within the associated region """
    if region['chr'] != chr:
        return True
    elif (start < region['start']) or (start > region['end']):
        return True
    elif (end <= region['start']) or (end > region['end']):
        return True
    return False


def load_coords(regions, coordFile=None, adjust=1, includeRegions=True):
    """ Reformat plot regions to required dictionary format and
        optionally load additional plot coordinates """
    coords = reformatRegions(regions)
    # Initalise empty dictionary with matching keys
    newCoords = {name : [] for name in coords.keys()}
    if coordFile is not None:
        with open(coordFile) as fh:
            for line in fh:
                chr, start, end, region = line.strip().split()
                start, end = adjustCoordinates(int(start), int(end), adjust)
                if region not in coords:
                    print(f'{region} not in regions file - skipping.',
                          file=sys.stderr)
                elif outOfRange(chr, start, end, regions.loc[region]):
                    print(f'Viewpoint {chr}:{start}-{end} not in '
                          f'range of region {region} - skipping.',
                          file=sys.stderr)
                else:
                    newCoords[region].append(f'{chr}_{start}_{end}')
    # If includeRegions is True then add original regions to coords
    if includeRegions:
        for region, pos in coords.items():
            newCoords[region].extend(pos)
    return newCoords
