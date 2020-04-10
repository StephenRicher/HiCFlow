#!/usr/bin/env python3

import re
import os
import sys

def main():
    for dir in sys.argv[1:]:
        # The lowest directory of dir is the sample name.
        sample = os.path.basename(dir.rstrip('/'))
        if '-' not in sample:
            sys.exit('Invalid sample name in path.')
        group = sample.split('-')[0]

        print(sample, dir, group, sep = '\t')

if __name__ == "__main__":
    main()
