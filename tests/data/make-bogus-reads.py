#!/usr/bin/env python

from khmer.utils import write_record
import screed
import sys

mutations = {
    0: (29, 19, 'T', 'A'),
    1: (19, 41, 'G', 'A'),
    2: (67, 5,  'G', 'T'),
    3: (63, 20, 'T', 'C'),
}

readnum = 0
for n, record in enumerate(screed.open(sys.argv[1])):
    if n in mutations:
        refrstart, mismatchpos, origbase, newbase = mutations[n]
        readseq = record.sequence[refrstart:refrstart+50]
        assert readseq[mismatchpos] == origbase
        mutseq = readseq[:mismatchpos] + newbase + readseq[mismatchpos+1:]
        oldseqname = record.name.split('-')[-1]

        readnum += 1
        record.name = 'read{}_{}_exact'.format(readnum, oldseqname)
        record.sequence = readseq
        record.quality = '3' * 50
        write_record(record, sys.stdout)

        readnum += 1
        record.name = 'read{}_{}_mismatch'.format(readnum, oldseqname)
        record.sequence = mutseq
        write_record(record, sys.stdout)
