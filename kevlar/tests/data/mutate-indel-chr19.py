#!/usr/bin/env python

from __future__ import print_function
from khmer.utils import write_record
import screed
import sys

mutations = {
    0: (42681, 10),
}

for n, record in enumerate(screed.open(sys.argv[1])):
    if n in mutations:
        start, dellength = mutations[n]
        seqlength = len(record.sequence)
        piece1 = record.sequence[:start]
        piece2 = record.sequence[start+dellength:]
        record.sequence = piece1 + piece2
        print('DEBUG ', piece1[-9:], '|', piece2[:9], sep='', file=sys.stderr)
        assert len(record.sequence) == seqlength - dellength

    write_record(record, sys.stdout)
