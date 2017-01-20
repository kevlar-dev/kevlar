#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict
import argparse
import sys

import khmer
from khmer.utils import write_record
from khmer import khmer_args
import kevlar
import screed


def subparser(subparsers):
    subparser = subparsers.add_parser('find')

    subparser.add_argument('--controls', metavar='FILE', nargs='+',
                           required=True, help='one or more countgraph files '
                           'corresponding to control sample(s)')
    subparser.add_argument('-x', '--ctrl_max', metavar='X', type=int,
                           default=1, help='k-mers with abund > X in any '
                           'control sample are uninteresting; default=1')
    subparser.add_argument('-y', '--case_min', metavar='Y',
                           type=int, default=5, help='ignore k-mers from case '
                           'with abund < Y; default=5')
    subparser.add_argument('-k', '--ksize', metavar='K', default=31, type=int,
                           help='k-mer size; default is 31')
    subparser.add_argument('-M', '--memory', default='1e6',
                           type=khmer_args.memory_setting, metavar='MEM',
                           help='total memory to allocate for each count '
                           'table; default is 1M')
    subparser.add_argument('--max-fpr', type=float, default=0.2, metavar='FPR',
                           help='terminate if the expected false positive rate'
                           ' is higher than the specified FPR; default is 0.2')
    subparser.add_argument('--out', type=argparse.FileType('w'),
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('--flush', action='store_true', help='flush output'
                           'after each read written')
    subparser.add_argument('--upint', type=float, default=1e6, metavar='INT',
                           help='debug update interval; default is 1000000')
    subparser.add_argument('--batch', type=int, nargs=2, metavar='INT',
                           help='process only a fraction of the k-mers in the '
                           'input file; invoke with "--batch X Y" where X is '
                           'a power of 2 and Y is between 1 and X inclusive')
    subparser.add_argument('case')


def load_case_and_controls(args):
    print('[kevlar::find] Loading case counttable', args.case, '...',
          end='', file=args.logfile)
    case = khmer.Counttable(args.ksize, args.memory / 4, 4)
    if args.batch:
        num_batches = int(args.batch[0])
        batch = int(args.batch[1])
        nr, nk = case.consume_fasta_banding(args.case, num_batches, batch - 1)
        message = 'batch {:d}/{:d} done'.format(batch, num_batches)
    else:
        message = 'done'
        nr, nk = case.consume_fasta(args.case)
    fpr = kevlar.calc_fpr(case)
    message += ', k={:d}'.format(case.ksize())
    message += '; {:d} reads and {:d} k-mers consumed'.format(nr, nk)
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    print(message, file=args.logfile)
    if fpr > args.max_fpr:
        sys.exit(1)

    controls = list()
    for ctlfile in args.controls:
        print('[kevlar::find] Loading control counttable', ctlfile, '...',
              end='', file=args.logfile)
        counttable = khmer.Counttable(args.ksize, args.memory / 4, 4)
        if args.batch:
            nr, nk = counttable.consume_fasta_banding(ctlfile, num_batches,
                                                      batch - 1)
            message = 'batch {:d}/{:d} done'.format(batch, num_batches)
        else:
            nr, nk = counttable.consume_fasta(ctlfile)
            message = 'done'
        assert counttable.ksize() == case.ksize()
        fpr = kevlar.calc_fpr(case)
        message += ', k={:d}'.format(counttable.ksize())
        message += '; {:d} reads and {:d} k-mers consumed'.format(nr, nk)
        message += '; estimated false positive rate is {:1.3f}'.format(fpr)
        print(message, file=args.logfile)
        if fpr > args.max_fpr:
            sys.exit(1)
        controls.append(counttable)

    return case, controls


def kmer_is_interesting(kmer, casecounts, controlcounts, case_min=5,
                        ctrl_max=1):
    caseabund = casecounts.get(kmer)
    if caseabund < case_min:
        return False

    ctrlabunds = list()
    for count in controlcounts:
        abund = count.get(kmer)
        if abund > ctrl_max:
            return False
        ctrlabunds.append(abund)

    return [caseabund] + ctrlabunds


def print_interesting_read(record, kmers, outstream, max_copy=2, flush=False):
    counts = defaultdict(int)
    for kmerpos in kmers:
        kmer = kmers[kmerpos][0]
        counts[kmer] += 1
        if counts[kmer] > max_copy:
            return

    write_record(record, outstream)
    for i in sorted(kmers):
        kmer = kmers[i][0]
        abunds = kmers[i][1:]
        abundstr = ' '.join([str(abund) for abund in abunds])
        print(' ' * i, kmer, ' ' * 10, abundstr, '#', sep='', file=outstream)
    if hasattr(outstream, 'flush') and flush:
        outstream.flush()


def main(args):
    case, controls = load_case_and_controls(args)

    print('[kevlar::find] Iterating over case reads', args.case,
          file=args.logfile)
    nkmers = 0
    nreads = 0
    unique_kmers = set()
    for n, record in enumerate(screed.open(args.case)):
        if n > 0 and n % args.upint == 0:
            print('    processed', n, 'reads...', file=args.logfile)

        read_novel_kmers = dict()
        for i, kmer in enumerate(case.get_kmers(record.sequence)):
            if args.batch:
                num_batches = int(args.batch[0])
                batch = int(args.batch[1]) - 1
                khash = case.hash(kmer)
                if khash & (num_batches - 1) != batch:
                    continue
            counts = kmer_is_interesting(kmer, case, controls, args.case_min,
                                         args.ctrl_max)
            if counts is False:
                continue
            read_novel_kmers[i] = [kmer] + counts
            minkmer = kevlar.revcommin(kmer)
            unique_kmers.add(minkmer)

        if len(read_novel_kmers) > 0:
            nreads += 1
            nkmers += len(read_novel_kmers)
            print_interesting_read(record, read_novel_kmers, args.out,
                                   args.flush)

    message = 'Found {:d} instances'.format(nkmers)
    message += ' of {:d} unique novel kmers'.format(len(unique_kmers))
    message += ' in {:d} reads'.format(nreads)
    print('[kevlar::find]', message, file=args.logfile)
