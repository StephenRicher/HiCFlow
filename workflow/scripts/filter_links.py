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
        '-p', '--p_value', default=0.05, type=float,
        help='P-value threshold for filtering interactions.')
    parser.add_argument(
        '-m', '--log_fc', default=0, type=float,
        help='Absolute log FC threshold for filtering interactions.')
    parser.add_argument(
        '--up', metavar='OUT_UP', help='Output file for UP LINKS file.')
    parser.add_argument(
        '--down', metavar='OUT_DOWN', help='Output file for DOWN LINKS file.')
    parser.set_defaults(function=filter_links)

    return (pct.execute(parser))

def filter_links(infile, up, down, log_fc, p_value):

    with ExitStack() as stack:

        f = stack.enter_context(pct.open(infile))
        out_up = stack.enter_context(pct.open(up, 'w'))
        out_down = stack.enter_context(pct.open(down, 'w'))

        for line in f:
            line_split = line.split()
            p_adj = float(line_split[9])
            m_adj = float(line_split[8])

            if abs(m_adj) >= log_fc and abs(p_adj) <= p_value:
                if m_adj > 0:
                    out_up.write(line)
                else:
                    out_down.write(line)

if __name__ == "__main__":
    main()
