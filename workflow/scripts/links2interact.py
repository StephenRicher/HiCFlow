#!/usr/bin/env python3

import sys
import pyCommonTools as pct
from contextlib import ExitStack

def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__, infile=True,
        in_type='LINKS')

    parser.add_argument(
        '--up', metavar='OUT_UP', help='Output file for UP interactions.')
    parser.add_argument(
        '--down', metavar='OUT_DOWN', help='Output file for DOWN interactions.')
    parser.set_defaults(function=links2interval)

    return (pct.execute(parser))

def links2interval(infile, up, down):

    sample1 = 'HB2_WT'
    sample2 = 'HB2_CL4'
    header = (f'track type=interact name="{sample1} vs {sample2} - all" '
              'description="Chromatin interactions - all" spectrum=on '
              'useScore=on maxHeightPixels=500:300:50 visibility=full')

    with ExitStack() as stack:

        f = stack.enter_context(pct.open(infile))
        out_up = stack.enter_context(pct.open(up, 'w'))
        out_down = stack.enter_context(pct.open(down, 'w'))

        if up == down:
            print(header,  file = out_up)
        else:
            print(header.replace('all', 'UP'), file = out_up)
            print(header.replace('all', 'DOWN'), file = out_down)

        for line in f:
            line = line.split()
            chr = line[0]
            start1 = line[1]
            end1 = line[2]
            start2 = line[4]
            end2 = line[5]
            # Retrive the M from index 8 that is NOT absolute value
            M = float(line[8])
            score = int(float(line[7]))

            if M > 0:
                out = out_up
            else:
                out = out_down

            print('chr1', start1, end2, '.', score, M, '.', 0,
                'chr1', start1, end1, '.', '.', 'chr1', start2, end2, '.', '.',
                sep = '\t', file = out)


if __name__ == "__main__":
    main()
