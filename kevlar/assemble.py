#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
import re
import sys


def assemble_fml_asm(partition, logstream=sys.stderr):
    reads = list(partition)
    assembler = kevlar.assembly.fml_asm(reads)
    for n, contig in enumerate(assembler, 1):
        name = 'contig{:d}'.format(n)
        record = kevlar.sequence.Record(name=name, sequence=contig)
        yield next(kevlar.augment.augment(reads, [record]))


def assemble(partstream, maxreads=10000, logstream=sys.stderr):
    n = 0
    pn = 0
    progress_indicator = kevlar.ProgressIndicator(
        '[kevlar::assemble] skipping partition with {counter} reads',
        interval=10, breaks=[100, 1000, 10000], usetimer=True,
        logstream=logstream,
    )
    for partid, partition in partstream:
        pn += 1
        progress_indicator.update()
        numreads = len(partition)
        if numreads > maxreads:  # pragma: no cover
            message = 'skipping partition with {:d} reads'.format(numreads)
            print('[kevlar::assemble] WARNING:', message, file=logstream)
            continue
        for contig in assemble_fml_asm(partition, logstream=logstream):
            n += 1
            newname = 'contig{}'.format(n)
            if partid is not None:
                newname += ' kvcc={}'.format(partid)
            contig.name = newname
            yield partid, contig
    message = 'processed {} partitions'.format(pn)
    message += ' and assembled {} contigs'.format(n)
    print('[kevlar::assemble]', message, file=logstream)


def main(args):
    readstream = kevlar.parse_augmented_fastx(kevlar.open(args.augfastq, 'r'))
    if args.part_id:
        pstream = kevlar.parse_single_partition(readstream, args.part_id)
    else:
        pstream = kevlar.parse_partitioned_reads(readstream)
    outstream = kevlar.open(args.out, 'w')
    assembler = assemble(pstream, maxreads=args.max_reads,
                         logstream=args.logfile)
    for partid, contig in assembler:
        kevlar.print_augmented_fastx(contig, outstream)
