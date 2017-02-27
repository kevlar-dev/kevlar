#!/usr/bin/env python

from __future__ import print_function
import sys

ahist = None
uhist = None

for infile in sys.argv[1:]:
    with open(infile, 'r') as fh:
        aline = next(fh)
        uline = next(fh)
        if ahist is None:
            ahist = eval(aline)
            uhist = eval(uline)
        else:
            ahist = list(map(sum, zip(ahist, eval(aline))))
            uhist = list(map(sum, zip(uhist, eval(uline))))
print(ahist)
print(uhist)
