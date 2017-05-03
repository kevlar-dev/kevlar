#!/usr/bin/env python

from __future__ import print_function
from khmer.utils import write_record
import screed
import sys

mutations = {
    0: (3566, 'A', 'C'),
    1: (9262, 'C', 'T'),
    2: (6134, 'C', 'G'),
}

for n, record in enumerate(screed.open(sys.argv[1])):
    if n in mutations:
        ind, origbase, newbase = mutations[n]
        assert record.sequence[ind] == origbase
        piece1 = record.sequence[:ind]
        piece2 = record.sequence[ind+1:]
        record.sequence = piece1 + newbase + piece2
        print('DEBUG ', piece1[-5:], '|', newbase, '|', piece2[:5], sep='',
              file=sys.stderr)

    write_record(record, sys.stdout)
