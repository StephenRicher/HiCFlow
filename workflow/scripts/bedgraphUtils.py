#!/usr/bin/env python3

import sys
import logging
import argparse
import pandas as pd
from collections import defaultdict


def splitScore(line, filetype):
    """ Split positions and score from BED/bedgraph record """

    assert filetype in ['bed', 'bedgraph']
    if filetype == 'bedgraph':
        chrom, start, end, score = line.split()
    else:
        try:
            chrom, start, end, name, score = line.split()[:5]
        except ValueError as e:
            logging.exception('Input BED file does not contain score column.')
    return chrom, int(start), int(end), float(score)


def readRegions(bed):
    """ Read BED file and return dict of chromosomes and allowed intervals """
    if bed is None:
        return None
    regions = defaultdict(list)
    with open(bed) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chrom, start, end = splitPos(line)
            regions[chrom].append(range(start, end))
    return regions


def hasOverlap(regions):
    """ Ensure there is no overlap in regions interval dictionary. """
    for chromIntervals in regions.values():
        allPos = []
        for interval in chromIntervals:
            allPos.extend(list(interval))
        if len(allPos) != len(set(allPos)):
            return True
    return False


def splitPos(line):
    chrom, start, end = line.split()[:3]
    return chrom, int(start), int(end)


def splitName(line):
    chrom, start, end, name = line.split()[:4]
    return chrom, int(start), int(end), name


class Bed:
    def __init__(self, line, buffer=0):
        self.record = line.strip('\n').split()
        self.buffer = int(buffer)
        assert len(self.record) >= 3

    @property
    def chrom(self):
        return self.record[0]

    @property
    def start(self):
        return int(self.record[1]) - self.buffer

    @property
    def end(self):
        return int(self.record[2]) + self.buffer

    @property
    def regionLength(self):
        return self.end - self.start

    @property
    def interval(self):
        return range(self.start, self.end)

    @property
    def name(self):
        try:
            return self.record[3]
        except IndexError:
            return None

    @property
    def score(self):
        try:
            return float(self.record[4])
        except IndexError:
            return None

    @property
    def strand(self):
        try:
            return self.record[5]
        except IndexError:
            return None

class Bedgraph:
    def __init__(self, line):
        self.record = line.strip('\n').split()
        assert len(self.record) >= 3

    @property
    def chrom(self):
        return self.record[0]

    @property
    def start(self):
        return int(self.record[1])

    @property
    def end(self):
        return int(self.record[2])

    @property
    def score(self):
        return float(self.record[3])

    @property
    def regionLength(self):
        return self.end - self.start

    @property
    def normScore(self):
        return self.score / self.regionLength

    @property
    def interval(self):
        return range(self.start, self.end)


def readBed(file, buffer=0, filetype='bed'):
    """ Construct sorted dictionary of Bed/Bedgraph objects """
    assert filetype in ['bed', 'bedgraph']
    records = defaultdict(list)
    with open(file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if filetype == 'bed':
                bed = Bed(line, buffer)
            else:
                bed = Bedgraph(line)
            records[bed.chrom].append(bed)
    # Sort per-chromosome ranges by start
    for chrom, bed in records.items():
        records[chrom] = sorted(bed, key=lambda r: r.start)
    return records


def readBedLength(file):
    """ Read BED as dictory of BED (keys) and BED length (values) """
    allEntries = {}
    with open(file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            entry = Bed(line)
            allEntries[entry] = entry.regionLength
    return allEntries
