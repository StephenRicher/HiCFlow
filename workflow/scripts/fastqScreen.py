#!/usr/bin/env python3

""" Wraper to run FastQScreen """

import sys
import glob
import shutil
import argparse
import tempfile
import subprocess
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'

def runFastQScreen(fastq, config, subset, threads, plotOut, dataOut):
    with tempfile.TemporaryDirectory() as dirname:
        command = ([
            'fastq_screen', '--outdir', dirname,
            '--force', '--aligner', 'bowtie2',
            '--conf', config, '--subset', str(subset),
            '--threads', str(threads), fastq
        ])
        subprocess.run(command)
        
        shutil.move(glob.glob(dirname + '/*_screen.txt')[0], dataOut)
        shutil.move(glob.glob(dirname + '/*_screen.png')[0], plotOut)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument('fastq', help='Input FastQ file.')
    parser.add_argument('config', help='Path to config file.')
    parser.add_argument(
        '--subset', type=int, default=100000,
        help='Number reads to sample (default: %(default)s)')
    parser.add_argument(
        '--threads', type=int, default=1,
        help='Threads for multiprocessing (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--plotOut', required=True, help='Path to save outplot plot.')
    requiredNamed.add_argument(
        '--dataOut', required=True, help='Path to save output text results.')
    parser.set_defaults(function=runFastQScreen)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
