#!/usr/bin/env python3

""" Process named-sorted SAM/BAM alignment files to identify fragment
    mapping, insert size and relative orientation of pairs.
"""

import re
import sys
import argparse
import fileinput
import pandas as pd
from typing import List
from collections import defaultdict
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def processHiC(sample: str, out: str, infile: List):


    data = defaultdict(list)

    with fileinput.input(infile) as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            try:
                read1 = Sam(line)
                read2 = Sam(next(fh))
            except StopIteration:
                return 1
            # Skip unmapped
            if (read1.rname == '*') or (read2.rname == '*'):
                continue
            read1, read2 = reorderReadPair(read1, read2)
            data['orientation'].append(getOrientation(read1, read2))
            data['cis'].append(isCisInteraction(read1, read2))
            data['readSeperation'].append(getReadSeperation(read1, read2))

    group, rep = sample.rsplit('-', 1)
    data = pd.DataFrame(data)
    data['group'] = group
    data['rep'] = rep
    data.to_pickle(out)


def getReadSeperation(read1, read2):
    return abs(read2.right_pos - read1.left_pos + 1)


def getOrientation(read1, read2):
    """
    Return relative orientation of read pairs. Assumes read pairs have
    been ordered such that read 1 is five prime of read 2.
    """

    if read1.is_reverse:
        if read2.is_reverse:
            orientation = "Same-reverse"
        else:
            orientation = "Outward"
    else:
        if read2.is_reverse:
            orientation = "Inward"
        else:
            orientation = "Same-forward"
    return orientation


def reorderReadPair(read1, read2):
    """
    Return a pair of reads such that read1 is left of read 2.
    Read pairs aligning to different chromosomes are returned unchanged.
    """

    if (isCisInteraction(read1, read2) and read1.left_pos > read2.left_pos):
        r1_reorder = read2
        r2_reorder = read1
    else:
        r1_reorder = read1
        r2_reorder = read2
    return r1_reorder, r2_reorder


def isCisInteraction(read1, read2):
    return read1.rname == read2.rname


class Sam:

    def __init__(self, record):
        record = record.strip().split('\t')
        self.flag = int(record[1])
        self.rname = record[2]
        self.left_pos = int(record[3])
        self.cigar = record[5]

    @property
    def is_reverse(self):
        return True if (self.flag & 0x10 != 0) else False

    @property
    def reference_length(self):
        cigar_split = re.findall(r'[A-Za-z]+|\d+', self.cigar)
        length = 0
        for idx, val in enumerate(cigar_split):
            if idx & 1 and val not in ["I", "S", "H", "P"]:
                length += int(cigar_split[idx-1])
        return length

    @property
    def right_pos(self):
        return self.left_pos + (self.reference_length - 1)


def parseArgs():

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=processHiC)
    parser.add_argument(
        'sample', help='Sample name in the format "group-rep".')
    parser.add_argument(
        'out', help='Path to write processed pickle file.')
    parser.add_argument(
        'infile', metavar='SAM', nargs='?', help='Input SAM file.')


    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
