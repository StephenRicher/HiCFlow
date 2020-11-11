#!/usr/bin/env python3

import argparse


def setDefaults(parser, verbose=True, version=None):
    """ Add version and verbose arguments to parser """

    if version:
        parser.add_argument('--version', action='version',
            version=f'%(prog)s {version}')
    if verbose:
        parser.add_argument(
            '--verbose', action='store_const', const=logging.DEBUG,
            default=logging.INFO, help='verbose logging for debugging')

        args = parser.parse_args()
        logFormat='%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
        logging.basicConfig(level=args.verbose, format=logFormat)
        del args.verbose
    else:
        args = parser.parse_args()

    return args
