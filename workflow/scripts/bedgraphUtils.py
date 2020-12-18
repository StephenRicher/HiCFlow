#!/usr/bin/env python3

import sys
import logging
import argparse
import pandas as pd
from collections import defaultdict


def splitScore(line, format):
    """ Split positions and score from BED/bedgraph record """

    assert format in ['bed', 'bedgraph']
    if format == 'bedgraph':
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
