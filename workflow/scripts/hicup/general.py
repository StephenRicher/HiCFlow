#!/usr/bin/env python3

""" Generic functions for running HiCUP wrapper scripts. """

__author__ = 'Stephen Richer'
__copyright__ = 'Copyright 2020, Stephen Richer'
__email__ = 'sr467@bath.ac.uk'
__license__ = 'MIT'
__version__ = '1.0.0'


import os
import re
import sys
import glob
import shutil
import pyCommonTools as pct


def get_filepath(dir, pattern):
    return glob.glob(os.path.join(dir, pattern))[0]


def move_file(source, destination=None, stderr=False):
    if destination is None:
        if stderr:
            destination = sys.stderr
        else:
            destination = sys.stdout
        with open(source) as f:
            shutil.copyfileobj(f, destination)
        os.remove(source)
    else:
        shutil.move(source, destination)


def set_zip(output, ext='.gz'):
    """ Return True is output file should be compressed. """

    log = pct.create_logger()
    zip_out = False
    if output is None:
        log.info('Uncompressed output will be written to stdout.')
    elif output.endswith(ext):
        zip_out = True
        log.info(f'Compressed output will be written to {output}.')
    else:
        log.info(f'Uncompressed output will be written to {output}.')

    return zip_out


def restriction_seq(value):
    """ Custom argument type for restriction enzyme argument. """

    if value.count('^') != 1:
        raise argparse.ArgumentTypeError(
            f'Restriction site {value} must contain one "^" at cut site.')
    elif re.search('[^ATCG^]', value, re.IGNORECASE):
        raise argparse.ArgumentTypeError(
            f'Restriction site {value} must only contain "ATCG^".')
    else:
        return value.upper()


def svg_file(value):
    """ Custom argument type for SVG file format. """

    if not value.endswith('.svg'):
        raise argparse.ArgumentTypeError(
            f'Output figure {value} must be end with ".svg"')
    else:
        return value
