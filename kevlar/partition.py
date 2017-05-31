#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import argparse
import sys
import networkx
import khmer
import kevlar
from kevlar.seqio import load_reads_and_kmers


def subparser(subparsers):
    subparser = subparsers.add_parser('partition')
    subparser.add_argument('-d', '--debug', action='store_true',
                           help='print debugging output')
    subparser.add_argument('-s', '--strict', action='store_true',
                           help='require perfect identity between overlapping '
                           'reads for inclusion in the same partition; by '
                           'default, only a shared interesting k-mer is '
                           'required')
    subparser.add_argument('--min-abund', metavar='X', type=int, default=2,
                           help='ignore k-mers with abundance lower than X; '
                           'default is 2')
    subparser.add_argument('--max-abund', metavar='Y', type=int, default=200,
                           help='ignore k-mers with abundance higher than Y; '
                           'default is 200')
    subparser.add_argument('--gml', metavar='FILE',
                           help='write read graph to .gml file')
    subparser.add_argument('outprefix', help='prefix for output files')
    subparser.add_argument('augfastq', help='original augmented Fastq file')


def main(args):
    debugout = None
    if args.debug:
        debugout = args.logfile

    timer = kevlar.Timer()
    timer.start()

    timer.start('loadreads')
    print('[kevlar::partition] Loading reads from', args.augfastq,
          file=args.logfile)
    inputfile = kevlar.open(args.augfastq, 'r')
    reads, kmers = load_reads_and_kmers(inputfile, debugout)
    elapsed = timer.stop('loadreads')
    print('[kevlar::partition]', 'Reads loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    timer.start('buildgraph')
    if args.strict:
        print('[kevlar::partition] Building read graph in strict mode',
              '(perfect match in read overlap required)', file=args.logfile)
        graph = kevlar.overlap.graph_init_strict(reads, kmers, args.min_abund,
                                                 args.max_abund, debugout)
    else:
        print('[kevlar::partition] Building read graph in relaxed mode',
              '(shared novel k-mer required)', file=args.logfile)
        graph = kevlar.overlap.graph_init_basic(kmers, logstream=debugout)
    elapsed = timer.stop('loadreads')
    print('[kevlar::partition]', 'Graph built in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    if args.gml:
        kevlar.to_gml(graph, args.gml, args.logfile)

    timer.start('writeoutput')
    print('[kevlar::partition] Writing output to prefix', args.outprefix,
          file=args.logfile)
    kevlar.overlap.write_partitions(graph, reads, args.outprefix, args.logfile)
    elapsed = timer.stop('writeoutput')
    print('[kevlar::partition]',
          'Output written in {:.2f} sec'.format(elapsed), file=args.logfile)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::partition]', message, file=args.logfile)
