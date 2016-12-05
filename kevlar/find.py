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

import khmer
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
    subparser.add_argument('--upint', type=float, default=1e6, help='debugging'
                           ' update interval; default is 1000000')
    subparser.add_argument('case_fastq')


def main(args):
    if args.kmers_out:
        unique_kmers = set()

    if args.paths_out:
        assert args.kmers_out, '--paths-out requires --kmers-out'
        linear_paths = defaultdict(set)

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

    print('[kevlar::find] Iterating over case reads', args.case_fastq,
          file=args.logfile)
    intreadcount = 0
    intkmercount = 0
    for n, record in enumerate(screed.open(args.case_fastq)):
        if n > 0 and n % args.upint == 0:
            print('    processed', n, 'reads...', file=args.logfile)
        novel_kmers = dict()
        for i, kmer in enumerate(case.get_kmers(record.sequence)):
            case_abund = case.get(kmer)
            if case_abund < args.case_lower_threshold:
                continue
            ctl_abunds = [c.get(kmer) for c in controls]
            ctl_thresh_pass = [a > args.control_threshold for a in ctl_abunds]
            if True in ctl_thresh_pass:
                continue

            novel_kmers[i] = [kmer, case_abund] + ctl_abunds
            intkmercount += 1
            if args.kmers_out:
                minkmer = revcommin(kmer)
                unique_kmers.add(minkmer)
                if args.paths_out:
                    linear_path = case.assemble_linear_path(minkmer)
                    min_path = revcommin(linear_path)
                    linear_paths[min_path].add(record.name)
        if len(novel_kmers) > 0:
            intreadcount += 1
            khmer.utils.write_record(record, args.out)
            for i in sorted(novel_kmers):
                kmer = novel_kmers[i][0]
                abunds = novel_kmers[i][1:]
                abundstr = ' '.join([str(abund) for abund in abunds])
                print(' ' * i, kmer, ' ' * 10, abundstr, '#', sep='',
                      file=args.out)
            if args.out == sys.stdout:
                sys.stdout.flush()

    if args.kmers_out:
        message = 'Found {} novel kmers in {} reads'.format(len(unique_kmers),
                                                            intreadcount)
        for i, kmer in enumerate(unique_kmers):
            kmername = 'kmer{}'.format(i+1)
            rcname = 'kmer{}_revcom'.format(i+1)
            rcseq = screed.dna.reverse_complement(kmer)
            print('>', kmer_name, '\n', kmer, sep='', file=args.kmers_out)
            print('>', rcname, '\n', rcseq, sep='', file=args.kmers_out)
        if args.paths_out:
            message += ', {} linear paths'.format(len(linear_paths))
            for i, linear_path in enumerate(linear_paths):
                readnames = lpaths[linear_path]
                lpathname = 'lpath{} {} reads {}'.format(
                                i+1, len(readnames),
                                ' '.join(readnames)
                            )
                rcname = 'lpath{}_revcom'.format(i+1)
                rcseq = screed.dna.reverse_complement(linear_path)
                print('>', lpathname, '\n', linear_path, sep='',
                      file=args.paths_out)
                print('>', rcname, '\n', rcseq, sep='', file=args.paths_out)
    else:
        message = 'Found {} novel kmers in {} reads'.format(intkmercount,
                                                            intreadcount)
    print('[kevlar::find]', message, file=args.logfile)
