#!/usr/bin/env python

import sys
import pyCommonTools as pct


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__, infile=True, in_type='GTF')
    parser.set_defaults(function=gff3_to_bed)
    return (pct.execute(parser))

def gff3_to_bed(infile):

    log = pct.create_logger()

    with pct.open(infile) as f:
        for index, line in enumerate(f):
            if line.startswith('#'):
                continue
            record = pct.GFF3(line)
            if record.feature != 'gene':
                continue
            try:
                name = record.attributes['Name']
            except KeyError:
                try:
                    name = record.attributes['ID'].strip('gene:')
                except KeyError:
                    log.error(f'No associated gene name or ID on line {index}')
                    continue

            print(
                record.seqname, record.start, record.end,
                name, 0, record.strand, sep='\t')

if __name__ == '__main__':
    main()
