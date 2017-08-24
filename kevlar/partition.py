#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import sys
import networkx
import khmer
import kevlar
from kevlar.seqio import load_reads_and_kmers


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
    graph = kevlar.ReadGraph()
    graph.load(inputfile, minabund=args.min_abund, maxabund=args.max_abund)
    elapsed = timer.stop('loadreads')
    print('[kevlar::partition]', 'Reads loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    timer.start('buildgraph')
    mode = 'strict' if args.strict else 'relaxed'
    message = 'Building read graph in {:s} mode'.format(mode)
    print('[kevlar::partition]', message, file=args.logfile)
    graph.populate_edges(strict=args.strict)
    elapsed = timer.stop('buildgraph')
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
