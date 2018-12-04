#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import networkx
import khmer
import kevlar


def partition(readstream, strict=False, minabund=None, maxabund=None,
              dedup=True, gmlfile=None):
    timer = kevlar.Timer()
    timer.start()

    timer.start('loadreads')
    kevlar.plog('[kevlar::partition] Loading reads')

    graph = kevlar.ReadGraph()
    graph.load(readstream, minabund=minabund, maxabund=maxabund)
    elapsed = timer.stop('loadreads')
    message = 'Reads loaded in {:.2f} sec'.format(elapsed)
    kevlar.plog('[kevlar::partition]', message)

    timer.start('buildgraph')
    mode = 'strict' if strict else 'relaxed'
    message = 'Building read graph in {:s} mode'.format(mode)
    kevlar.plog('[kevlar::partition]', message)
    graph.populate_edges(strict=strict)
    elapsed = timer.stop('buildgraph')
    message = 'Graph built in {:.2f} sec'.format(elapsed)
    kevlar.plog('[kevlar::partition]', message)

    if gmlfile:  # pragma: no cover
        kevlar.to_gml(graph, gmlfile, logstream)

    timer.start('partition')
    kevlar.plog('[kevlar::partition] Partition readgraph')
    part_iter = graph.partitions(dedup, minabund, maxabund, abundfilt=True)
    for n, part in enumerate(part_iter, 1):
        reads = [graph.get_record(readname) for readname in list(part)]
        for read in reads:
            read.name += ' kvcc={:d}'.format(n)
        yield n, reads
    elapsed = timer.stop('partition')
    message = 'Partitioning done in {:.2f} sec'.format(elapsed)
    kevlar.plog('[kevlar::partition]', message)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    kevlar.plog('[kevlar::partition]', message)


def main(args):
    if args.split:
        kevlar.mkdirp(args.split, trim=True)
    outstream = None if args.split else kevlar.open(args.out, 'w')
    readstream = kevlar.parse_augmented_fastx(kevlar.open(args.infile, 'r'))
    partitioner = partition(
        readstream, strict=args.strict, minabund=args.min_abund,
        maxabund=args.max_abund, dedup=args.dedup, gmlfile=args.gml,
    )
    numreads = 0
    for partnum, part in partitioner:
        numreads += len(part)
        if args.split:
            ofname = '{:s}.cc{:d}.augfastq.gz'.format(args.split, partnum)
            with kevlar.open(ofname, 'w') as outfile:
                for read in part:
                    kevlar.print_augmented_fastx(read, outfile)
        else:
            for read in part:
                kevlar.print_augmented_fastx(read, outstream)
    message = 'grouped {:d} reads'.format(numreads)
    message += ' into {:d} connected components'.format(partnum)
    kevlar.plog('[kevlar::partition]', message)
