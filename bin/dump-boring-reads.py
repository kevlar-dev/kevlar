#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys

from Bio import SeqIO
import pysam
import khmer


def match(record, seq):
    refrseq = str(seq[record.pos:record.pos+record.rlen].seq)
    return refrseq == records.seq


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seqid', metavar='SEQ',
                        help='dump reads not mapped to SEQ')
    parser.add_argument('--genomemask', metavar='FILE', help='dump reads with '
                        'median k-mer abundance >= 1 in the specified genome; '
                        'if both --seqid and --genomemask are declared, reads '
                        'passing either filter will be kept')
    parser.add_argument('--maskmemory', metavar='SIZE', default=2e9, type=float,
                        help='memory to be occupied by genome mask; default is'
                        ' 2e9 (2G)')
    parser.add_argument('--mask-k', metavar='K', default=31, type=int,
                        help='k size for genome mask')
    parser.add_argument('--out', type=argparse.FileType('w'),
                        help='output file; default is terminal (stdout)')
    parser.add_argument('--dry-run', action='store_true',
                        help='do not execute, just print arguments')
    parser.add_argument('fasta')
    parser.add_argument('bam')
    return parser


def main(args):
    if args.dry_run:
        print(*sys.argv)
        exit(0)

    with open(args.fasta, 'r') as seqfile:
        seqs = {record.id: record for record in SeqIO.parse(seqfile, 'fasta')}

    if args.genomemask:
        genomemask = khmer.Countgraph(args.mask_k, int(args.maskmemory / 4), 4)
        genomemask.consume_fasta(args.genomemask)

    bam = pysam.AlignmentFile(args.bam, 'rb')
    for i, record in enumerate(bam):
        if (i+1) % 10000 == 0:
            print('...processed', i, 'records', file=sys.stderr)

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
        if not pass_seqid_filter and not pass_genomemask_filter:
            continue

        # Discard reads that match the reference genome exactly
        matchcigar = '{:d}M'.format(record.rlen)
        if record.cigarstring == matchcigar:
            seq = seqs[bam.get_reference_name(record.tid)]
            refrseq = str(seq[record.pos:record.pos+record.rlen].seq)
            if refrseq.upper() == record.seq.upper():
                continue

        print('@', record.qname, '\n', record.seq, '\n+\n', record.qual,
              sep='', file=args.out)

if __name__ == '__main__':
    main(get_parser().parse_args())
