#!/usr/bin/env python

import sys
import straw

result = straw.straw(
    sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6]))

out = sys.argv[7]
with open(out, 'w') as f:
    for i in range(len(result[0])):
        print(result[0][i], result[1][i], result[2][i], sep = '\t', file = f)

