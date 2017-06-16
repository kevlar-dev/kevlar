#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import pysam
import khmer
from khmer import khmer_args

import kevlar


def match(record, seq):
    refrseq = str(seq[record.pos:record.pos+record.rlen].seq)
    return refrseq == records.seq


def main(args):
    print('[kevlar::dump] Loading reference sequence', file=args.logfile)
    with kevlar.open(args.refr, 'r') as genome:
        seqs = kevlar.seqio.parse_seq_dict(genome)

    if args.genomemask:
        print('[kevlar::dump] Loading genome mask', file=args.logfile)
        genomemask = khmer.Countgraph(args.mask_k, int(args.maskmemory / 4), 4)
        genomemask.consume_seqfile(args.genomemask)

    bam = pysam.AlignmentFile(args.reads, 'rb')
    fastq = kevlar.open(args.out, 'w')
    for i, record in enumerate(bam):
        if i > 0 and i % 50000 == 0:
            print('...processed', i, 'records', file=args.logfile)

        # We only want primary sequences
        if record.is_secondary or record.is_supplementary:
            continue

        # If we're restricting to a particular chromosome, keep reads with that
        # chromosome's seqid or reads whose median k-mer abundance is 0 in the
        # genome mask.
        pass_seqid_filter = True
        if args.seqid:
            if record.reference_id < 0 or \
               bam.get_reference_name(record.reference_id) != args.seqid:
                pass_seqid_filter = False
        pass_genomemask_filter = True
        if args.genomemask:
            count = genomemask.get_median_count(record.seq)[0]
            pass_genomemask_filter = count < 1

        if args.seqid and args.genomemask:
            if not pass_seqid_filter and not pass_genomemask_filter:
                continue
        elif args.seqid:
            if not pass_seqid_filter:
                continue
        elif args.genomemask:
            if not pass_genomemask_filter:
                continue

        # Discard reads that match the reference genome exactly
        matchcigar = '{:d}M'.format(record.rlen)
        if record.cigarstring == matchcigar:
            seq = seqs[bam.get_reference_name(record.tid)]
            refrseq = str(seq[record.pos:record.pos+record.rlen])
            if refrseq.upper() == record.seq.upper():
                continue

        # Make sure paired reads have unique IDs
        name = record.qname
        if record.flag & 1:
            # Logical XOR: if the read is paired, it should be first in pair
            # or second in pair, not both.
            assert (record.flag & 64) != (record.flag & 128)
            suffix = '/1' if record.flag & 64 else '/2'
            if not name.endswith(suffix):
                name += suffix

        print('@', name, '\n', record.seq, '\n+\n', record.qual, sep='',
              file=fastq)
