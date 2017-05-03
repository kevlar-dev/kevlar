#!/usr/bin/env python

from __future__ import print_function
from khmer.utils import write_record
import screed
import sys

mutations = {
    0: (11015, 5),
}

for n, record in enumerate(screed.open(sys.argv[1])):
    if n in mutations:
        start, dellength = mutations[n]
        seqlength = len(record.sequence)
        piece1 = record.sequence[:start]
        piece2 = record.sequence[start+dellength:]
        record.sequence = piece1 + piece2
        print('DEBUG ', piece1[-5:], '|', piece2[:5], sep='', file=sys.stderr)
        assert len(record.sequence) == seqlength - dellength

    write_record(record, sys.stdout)
