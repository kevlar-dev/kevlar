#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from itertools import cycle
import kevlar


def split(pstream, outstreams, maxreads=10000):
    """Split the partitions across the N outstreams."""
    progress_indicator = kevlar.ProgressIndicator(
        '[kevlar::split] processed {counter} partitions',
        interval=100, breaks=[1000, 10000, 100000], usetimer=True,
    )
    for partdata, outstream in zip(pstream, cycle(outstreams)):
        partid, partition = partdata
        if len(partition) > maxreads:
            message = 'WARNING: discarding partition '
            message += 'with {} reads'.format(len(partition))
            kevlar.plog('[kevlar::split]', message)
            continue
        for read in partition:
            kevlar.print_augmented_fastx(read, outstream)
        progress_indicator.update()


def main(args):
    partfile = kevlar.open(args.infile, 'r')
    readstream = kevlar.parse_augmented_fastx(partfile)
    partstream = kevlar.parse_partitioned_reads(readstream)
    outstreams = list()
    for i in range(args.numfiles):
        outfile = '{:s}.{:d}.augfastx'.format(args.base, i)
        if args.infile.endswith('.gz'):
            outfile += '.gz'
        os = kevlar.open(outfile, 'w')
        outstreams.append(os)
    split(partstream, outstreams)
