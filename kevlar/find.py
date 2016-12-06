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
from sys import stdout

import khmer
from khmer.utils import write_record
import kevlar
import screed


def subparser(subparsers):
    subparser = subparsers.add_parser('find')

    subparser.add_argument('--controls', metavar='FILE', nargs='+',
                           required=True, help='one or more countgraph files '
                           'corresponding to control sample(s)')
    subparser.add_argument('--case', metavar='FILE', required=True,
                           help='countgraph of sample from case/condition of '
                           'interest')
    subparser.add_argument('-x', '--ctrl_max', metavar='X', type=int,
                           default=1, help='k-mers with abund > X in any '
                           'control sample are uninteresting; default=1')
    subparser.add_argument('-y', '--case_min', metavar='Y',
                           type=int, default=5, help='ignore k-mers from case '
                           'with abund < Y; default=5')
    subparser.add_argument('--out', type=argparse.FileType('w'),
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('--kmers-out', type=argparse.FileType('w'),
                           default=None, metavar='FILE',
                           help='output novel k-mers to specified file')
    subparser.add_argument('--paths-out', type=argparse.FileType('w'),
                           default=None, metavar='FILE',
                           help='output linear paths to specified file')
    subparser.add_argument('--upint', type=float, default=1e6, metavar='INT',
                           help='debug update interval; default is 1000000')
    subparser.add_argument('case_fastq')


def load_case_and_controls(args):
    print('[kevlar::find] Loading case countgraph', args.case, '...',
          end='', file=args.logfile)
    case = khmer.load_countgraph(args.case)
    print('done, k={:d}'.format(case.ksize()), file=args.logfile)

    controls = list()
    for ctlfile in args.controls:
        print('[kevlar::find] Loading control countgraph', ctlfile, '...',
              end='', file=args.logfile)
        countgraph = khmer.load_countgraph(ctlfile)
        assert countgraph.ksize() == case.ksize()
        controls.append(countgraph)
        print('done', file=args.logfile)

    return case, controls


def main(args):
    case, controls = load_case_and_controls(args)

    variants = kevlar.VariantSet()
    print('[kevlar::find] Iterating over case reads', args.case_fastq,
          file=args.logfile)
    for n, record in enumerate(screed.open(args.case_fastq)):
        if n > 0 and n % args.upint == 0:
            print('    processed', n, 'reads...', file=args.logfile)
        read_novel_kmers = dict()
        for i, kmer in enumerate(case.get_kmers(record.sequence)):
            case_abund = case.get(kmer)
            if case_abund < args.case_min:
                continue
            ctl_abunds = [c.get(kmer) for c in controls]
            ctl_thresh_pass = [a > args.ctrl_max for a in ctl_abunds]
            if True in ctl_thresh_pass:
                continue

            lpath = case.assemble_linear_path(kmer)
            variants.add_kmer(kmer, record.name, lpath)
            read_novel_kmers[i] = [kmer, case_abund] + ctl_abunds

        if len(read_novel_kmers) > 0:
            write_record(record, args.out)
            for i in sorted(read_novel_kmers):
                kmer = read_novel_kmers[i][0]
                abunds = read_novel_kmers[i][1:]
                abundstr = ' '.join([str(abund) for abund in abunds])
                print(' ' * i, kmer, ' ' * 10, abundstr, '#', sep='',
                      file=args.out)
            if args.out == stdout:
                stdout.flush()

    if args.kmers_out:
        variants.kmer_table(args.kmers_out)
    if args.paths_out:
        variants.path_table(args.paths_out)

    message = 'Found {} novel kmers in {} reads, {} linear paths'.format(
        variants.nkmers, variants.nreads, variants.npaths
    )
    print('[kevlar::find]', message, file=args.logfile)
