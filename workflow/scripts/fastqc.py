#!/usr/bin/env python3

""" Wraper to run FastQC """

import sys
import glob
import shutil
import argparse
import tempfile
import subprocess
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'

def runFastQC(fastq, htmlOut, dataOut):
    with tempfile.TemporaryDirectory() as dirname:
        command = ['fastqc', '--outdir', dirname, fastq]
        subprocess.run(command)
        shutil.move(glob.glob(dirname + '/*_fastqc.zip')[0], dataOut)
        shutil.move(glob.glob(dirname + '/*_fastqc.html')[0], htmlOut)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument('fastq', help='Input FastQ file.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--htmlOut', required=True, help='Path to save outplot html.')
    requiredNamed.add_argument(
        '--dataOut', required=True, help='Path to save output text results.')
    parser.set_defaults(function=runFastQC)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
