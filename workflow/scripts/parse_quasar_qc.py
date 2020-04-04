#!/usr/bin/env python3

import sys
import pyCommonTools as pct

def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__, infile=True,
        in_type='QUASAR_QC', nargs='*')
    parser.set_defaults(function=parse_quasar)

    return (pct.execute(parser))

def parse_quasar(infile):
    for n, report in enumerate(infile):
        with pct.open(report) as f:
            in_qc_block = False
            for line in f:
                line = line.strip()
                if line.startswith('Resolution'):
                    in_qc_block = True
                    ref = line.split()[3]
                    if n == 0:
                        sys.stdout.write(
                            'resolution\tcoverage\tall_qc\tref_qc\tref\n')
                elif in_qc_block:
                    if not line:
                        break
                    else:
                        line = line.split()
                        res = int(line[0]) * 1000
                        res = res * 1000 if line[1].lower() == 'mb' else res
                        cov = line[2].replace(',','')
                        all_qc = line[3]
                        ref_qc = line[4]
                        sys.stdout.write(
                            f'{res}\t{cov}\t{all_qc}\t{ref_qc}\t{ref}\n')


if __name__ == "__main__":
    main()
