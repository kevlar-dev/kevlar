#!/usr/bin/env python

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
        record.sequence = record.sequence[:start] + record.sequence[start+dellength:]
        assert len(record.sequence) == seqlength - dellength

    write_record(record, sys.stdout)
