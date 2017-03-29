#!/usr/bin/env python

from khmer.utils import write_record
import screed
import sys

mutations = {
    0: [(65, 21, 10)],
    1: [(304, 5, 30), (3025, 24, 150)],
    2: [(964, 43, 120)],
}

readnum = 0
for n, record in enumerate(screed.open(sys.argv[1])):
    if n in mutations:
        origseq = record.sequence
        oldseqname = record.name.split('-')[-1]
        for refrpos, indelpos, dellength in mutations[n]:
            length1 = indelpos
            length2 = 50 - length1
            start2 = refrpos + indelpos + dellength
            part1 = record.sequence[refrpos:refrpos+indelpos]
            part2 = record.sequence[start2:start2+length2]

            readnum += 1
            readname = 'read{}_{}_{}bpindel'.format(readnum, oldseqname,
                                                    dellength)
            record.sequence = part1 + part2
            record.name = readname
            assert len(record.sequence) == 50
            record.quality = '3' * 50
            write_record(record, sys.stdout)
            record.sequence = origseq
            del record.quality
