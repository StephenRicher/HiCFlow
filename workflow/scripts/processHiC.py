#!/usr/bin/env python3

""" Process named-sorted SAM/BAM alignment files to identify fragment
    mapping, insert size, ditag size and relative orientation of pairs.
"""

import os
import re
import sys
import bisect
import argparse
import fileinput
from typing import List
from collections import defaultdict
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def processHiC(infile: List, digest: str):

    print('sample', 'orientation', 'interaction_type', 'ditag_length', 'insert_size', sep='\t')

    digest = processDigest(digest)
    refSizes = defaultdict(int)

    with fileinput.input(infile) as fh:
        for line in fh:
            if line.startswith('@'):
                if line.startswith('@SQ'):
                    ref, size = line.split()[1:]
                    ref = ref.split(':')[1]
                    size = int(size.split(':')[1])
                    refSizes[ref] = size
            else:
                try:
                    read1 = Sam(line)
                    read2 = Sam(next(fh))
                except StopIteration:
                    return 1
                # Skip unmapped
                if (read1.rname not in digest) or (read2.rname not in digest):
                    continue
                read1, read2 = reorderReadPair(read1, read2)
                print(trimPath(infile), end='\t')
                print(getOrientation(read1, read2), end='\t')
                print(interactionType(read1, read2), end='\t')
                print(getDitagLength(read1, read1, digest, refSizes), end='\t')
                print(getInsertSize(read1, read2), end='\n')


def trimPath(path, suffix='.subsampled.sam'):
    """ Remove file path and defined suffix """
    base = os.path.basename(path)
    return re.sub(f'{suffix}$', '', base)


def getDitagLength(read1, read2, digest, refSizes):
    read1Frag = getFragment(read1, digest, refSizes)
    read2Frag = getFragment(read2, digest, refSizes)
    read1TagLength = tagLength(read1, read1Frag)
    read2TagLength = tagLength(read2, read2Frag)
    return read1TagLength + read2TagLength


def getInsertSize(read1, read2):
    return read2.right_pos - read1.left_pos + 1


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

    if (interactionType(read1, read2) == "cis"
            and read1.left_pos > read2.left_pos):
        r1_reorder = read2
        r2_reorder = read1
    else:
        r1_reorder = read1
        r2_reorder = read2
    return r1_reorder, r2_reorder


def getFragment(read, digest, refSizes):
    frag = bisect.bisect_left(digest[read.rname], read.middle_pos)
    start = 1 if frag == 0 else digest[read.rname][frag - 1] + 1
    try:
        end = digest[read.rname][frag]
    except IndexError:
        end = refSizes[read.rname]
    return fragment(frag, start, end)


def interactionType(read1, read2):
    if read1.rname != read2.rname:
        interaction = "trans"
    else:
        interaction = "cis"
    return interaction


def tagLength(read, fragment):
    if read.is_reverse:
        return read.five_prime_pos - fragment.start + 1
    else:
        return fragment.end - read.five_prime_pos + 1


class fragment:
    def __init__(self, number, start, end):
        self.number = number
        self.start = start
        self.end = end


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
    def is_read1(self):
        return True if (self.flag & 0x40 != 0) else False

    @property
    def is_paired(self):
        return True if (self.flag & 0x1 != 0) else False

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

    @property
    def five_prime_pos(self):
        if self.is_reverse:
            return self.right_pos
        else:
            return self.left_pos

    @property
    def three_prime_pos(self):
        return self.left_pos if self.is_reverse else self.right_pos

    @property
    def middle_pos(self):
        return round((self.right_pos + self.left_pos)/2)


def processDigest(digest):
    d = defaultdict(list)
    with open(digest) as fh:
        for fragment in fh:
            ref, start = fragment.split()[:2]
            d[ref].append(int(start))
    return(d)


def parseArgs():

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=processHiC)
    parser.add_argument(
        'digest', metavar='DIGEST',
        help='Digest file from HiCExplorer findRestSites')
    parser.add_argument(
        'infile', metavar='SAM', nargs='?', help='Input SAM file.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
