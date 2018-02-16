#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import json
import sys
import networkx
import khmer
import kevlar
import screed


def load_mate_map(matefile, mapfile=None):
    """Load mate sequences into a dictionary.

    The `matefile` variable contains mates of interesting reads in Fast[aq]
    format, and `mapfile` contains a mapping of interesting read IDs to mate
    read IDs in JSON. Returns a dictionary with interesting read IDs as keys
    and the corresponding mate as a read record.
    """
    if mapfile is None:
        mapfile = matefile + '.map.json'
    matemap = json.load(open(mapfile, 'r'))
    matereads = dict()
    for record in screed.open(matefile):
        matereads[record.name] = record
    return matereads, matemap


def partition(readstream, strict=False, minabund=None, maxabund=None,
              dedup=True, gmlfile=None, logstream=sys.stderr):
    timer = kevlar.Timer()
    timer.start()

    timer.start('loadreads')
    print('[kevlar::partition] Loading reads', file=logstream)

    graph = kevlar.ReadGraph()
    graph.load(readstream, minabund=minabund, maxabund=maxabund)
    elapsed = timer.stop('loadreads')
    print('[kevlar::partition]', 'Reads loaded in {:.2f} sec'.format(elapsed),
          file=logstream)

    timer.start('buildgraph')
    mode = 'strict' if strict else 'relaxed'
    message = 'Building read graph in {:s} mode'.format(mode)
    print('[kevlar::partition]', message, file=logstream)
    graph.populate_edges(strict=strict)
    elapsed = timer.stop('buildgraph')
    print('[kevlar::partition]', 'Graph built in {:.2f} sec'.format(elapsed),
          file=logstream)

    if gmlfile:
        kevlar.to_gml(graph, gmlfile, logstream)

    timer.start('partition')
    print('[kevlar::partition] Partition readgraph', file=logstream)
    part_iter = graph.partitions(dedup, minabund, maxabund, abundfilt=True)
    for part in part_iter:
        reads = [graph.get_record(readname) for readname in list(part)]
        yield reads
    elapsed = timer.stop('partition')
    print('[kevlar::partition]',
          'Partitioning done in {:.2f} sec'.format(elapsed), file=logstream)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::partition]', message, file=logstream)


def main(args):
    if args.split:
        kevlar.mkdirp(args.split, trim=True)
    matereads, matemap, mateout = None, None, None
    if args.mate_file:
        if not args.split and not args.mate_out:
            msg = 'must declare "--split" or "--mate-out" with "--mate-file"'
            raise ValueError(msg)
        if args.mate_out:
            if args.split:
                msg = 'WARNING: ignoring mate output file in "--split" mode'
                print('[kevlar::partition]', msg, file=args.logfile)
            else:
                mateout = kevlar.open(args.mate_out, 'w')
        matereads, matemap = load_mate_map(args.mate_file)
    outstream = None if args.split else kevlar.open(args.out, 'w')
    readstream = kevlar.parse_augmented_fastx(kevlar.open(args.infile, 'r'))
    partitioner = partition(readstream, strict=args.strict,
                            minabund=args.min_abund, maxabund=args.max_abund,
                            dedup=args.dedup, gmlfile=args.gml,
                            logstream=args.logfile)
    partnum = 0
    numreads = 0
    for partnum, part in enumerate(partitioner, 1):
        numreads += len(part)
        if args.split:
            ofname = '{:s}.cc{:d}.augfastq.gz'.format(args.split, partnum)
            if args.mate_file:
                mateoutfile = '{:s}.cc{:d}.mates.augfastq.gz'.format(
                    args.split, partnum
                )
                mateout = kevlar.open(mateoutfile, 'w')
            with kevlar.open(ofname, 'w') as outfile:
                for read in part:
                    kevlar.print_augmented_fastx(read, outfile)
                    if args.mate_file:
                        mateid = matemap[read.name]
                        materead = matereads[mateid]
                        khmer.utils.write_record(materead, mateout)

        else:
            for read in part:
                if mateout:
                    mateid = matemap[read.name]
                    materead = matereads[mateid]
                    materead.name += ' kvcc={:d}'.format(partnum)
                    khmer.utils.write_record(materead, mateout)
                read.name += ' kvcc={:d}'.format(partnum)
                kevlar.print_augmented_fastx(read, outstream)
    message = '[kevlar::partition] grouped {:d} reads'.format(numreads)
    message += ' into {:d} connected components'.format(partnum)
    print(message, file=args.logfile)
