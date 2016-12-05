#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import argparse

import pysam
import khmer

import kevlar


def match(record, seq):
    refrseq = str(seq[record.pos:record.pos+record.rlen].seq)
    return refrseq == records.seq


def subparser(subparsers):
    subparser = subparsers.add_parser('dump')
    subparser.add_argument('--seqid', metavar='SEQ',
                           help='dump reads not mapped to SEQ')
    subparser.add_argument('--genomemask', metavar='FILE', help='dump reads '
                           'with median k-mer abundance >= 1 in the specified '
                           'genome; if both --seqid and --genomemask are '
                           'declared, reads passing either filter will be '
                           'kept')
    subparser.add_argument('--maskmemory', metavar='SIZE', default=2e9,
                           type=float, help='memory to be occupied by genome '
                           'mask; default is 2e9 (2G)')
    subparser.add_argument('--mask-k', metavar='K', default=31, type=int,
                           help='k size for genome mask')
    subparser.add_argument('--out', metavar='FILE',
                           type=argparse.FileType('w'),
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('refr', help='reference sequence in Fasta format')
    subparser.add_argument('reads', help='read alignments in BAM format')


def main(args):
    print('[kevlar::dump] Loading reference sequence', file=args.logfile)
    with open(args.refr, 'r') as genome:
        seqs = kevlar.fasta.parse_seq_dict(genome)

    if args.genomemask:
        print('[kevlar::dump] Loading genome mask', file=args.logfile)
        genomemask = khmer.Countgraph(args.mask_k, int(args.maskmemory / 4), 4)
        genomemask.consume_fasta(args.genomemask)

    bam = pysam.AlignmentFile(args.reads, 'rb')
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

        print('@', record.qname, '\n', record.seq, '\n+\n', record.qual,
              sep='', file=args.out)
