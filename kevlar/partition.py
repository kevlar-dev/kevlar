#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import networkx
import khmer
import kevlar
from kevlar.seqio import load_reads_and_kmers


def write_partitions(graph, ccprefix, dedup=True, minabund=None, maxabund=None,
                     logstream=sys.stderr):
    """Given a read graph, write distinct partitions to separate files."""
    n = 0
    reads_in_ccs = 0
    cclog = open(ccprefix + '.cc.log', 'w')
    part_iter = graph.partitions(dedup, minabund, maxabund)
    for n, cc in enumerate(part_iter):
        readnames = [r for r in cc]
        print('CC', n, len(cc), readnames, sep='\t', file=cclog)
        reads_in_ccs += len(cc)
        outfilename = '{:s}.cc{:d}.augfastq.gz'.format(ccprefix, n)
        with kevlar.open(outfilename, 'w') as outfile:
            for readid in cc:
                record = graph.get_record(readid)
                kevlar.print_augmented_fastx(record, outfile)
    message = '[kevlar::partition] grouped {:d} reads'.format(reads_in_ccs)
    message += ' into {:d} connected components'.format(n + 1)
    print(message, file=logstream)


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
    readstream = kevlar.parse_augmented_fastx(inputfile)
    graph = kevlar.ReadGraph()
    graph.load(readstream, minabund=args.min_abund, maxabund=args.max_abund)
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
    write_partitions(graph, args.outprefix, dedup=args.dedup,
                     minabund=args.min_abund, maxabund=args.max_abund,
                     logstream=args.logfile)
    elapsed = timer.stop('writeoutput')
    print('[kevlar::partition]',
          'Output written in {:.2f} sec'.format(elapsed), file=args.logfile)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::partition]', message, file=args.logfile)
