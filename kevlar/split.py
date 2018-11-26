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


def split(pstream, outstreams):
    """Split the partitions across the N outstreams."""
    for partdata, outstream in zip(pstream, cycle(outstreams)):
        partid, partition = partdata
        for read in partition:
            kevlar.print_augmented_fastx(read, outstream)


def main(args):
    partfile = kevlar.open(args.infile, 'r')
    readstream = kevlar.parse_augmented_fastx(partfile)
    partstream = kevlar.parse_partitioned_reads(readstream)
    outstreams = list()
    for i in range(args.numfiles):
        outfile = '{:s}.{:d}'.format(args.base, i + 1)
        os = kevlar.open(outfile, 'w')
        outstreams.append(os)
    split(partstream, outstreams)
