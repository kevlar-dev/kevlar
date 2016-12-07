#!/usr/bin/env python

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
        record.sequence = record.sequence[:ind] + newbase + record.sequence[ind+1:]

    write_record(record, sys.stdout)
