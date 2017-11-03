#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import pysam
import screed
import kevlar
from khmer.utils import write_record


def perfectmatch(record, bam, refrseqs):
    matchcigar = '{:d}M'.format(record.rlen)
    if record.cigarstring == matchcigar:
        seq = refrseqs[bam.get_reference_name(record.tid)]
        refrsubseq = str(seq[record.pos:record.pos+record.rlen])
        if refrsubseq.upper() == record.seq.upper():
            return True
    return False


def readname(record):
    name = record.qname
    if record.flag & 1:
        # Logical XOR: if the read is paired, it should be first in pair
        # or second in pair, not both.
        assert (record.flag & 64) != (record.flag & 128)
        suffix = '/1' if record.flag & 64 else '/2'
        if not name.endswith(suffix):
            name += suffix
    return name


def dump(bamstream, refrstream, upint=50000, logstream=sys.stderr):
    print('[kevlar::dump] Loading reference sequence', file=logstream)
    refrseqs = kevlar.seqio.parse_seq_dict(refrstream)

    bam = pysam.AlignmentFile(bamstream, 'rb')
    for i, record in enumerate(bam, 1):
        if i % upint == 0:
            print('...processed', i, 'records', file=logstream)
        if record.is_secondary or record.is_supplementary:
            continue
        if perfectmatch(record, bam, refrseqs):
            continue
        rn = readname(record)
        yield screed.Record(name=rn, sequence=record.seq, quality=record.qual)


def main(args):
    fastq = kevlar.open(args.out, 'w')
    refrstream = kevlar.open(args.refr, 'r')
    for read in dump(args.reads, refrstream, logstream=args.logfile):
        write_record(read, fastq)
