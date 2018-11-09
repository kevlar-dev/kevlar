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
    ccmatch = re.search(r'kvcc=(\d+)', reads[0].name)
    cc = ccmatch.group(1) if ccmatch else None
    assembler = kevlar.assembly.fml_asm(reads)
    for n, contig in enumerate(assembler, 1):
        name = 'contig{:d}'.format(n)
        if cc is not None:
            name += ' kvcc={}'.format(cc)
        record = kevlar.sequence.Record(name=name, sequence=contig)
        yield next(kevlar.augment.augment(reads, [record]))


def assemble(partstream, logstream=sys.stderr):
    n = 0
    for partition in partstream:
        for contig in assemble_fml_asm(partition, logstream=logstream):
            n += 1
            contig.name = re.sub(r'^contig\d+', 'contig' + str(n), contig.name)
            yield contig


def main(args):
    readstream = kevlar.parse_augmented_fastx(kevlar.open(args.augfastq, 'r'))
    if args.part_id:
        pstream = kevlar.parse_single_partition(readstream, args.part_id)
    else:
        pstream = kevlar.parse_partitioned_reads(readstream)
    outstream = kevlar.open(args.out, 'w')
    for contig in assemble(pstream, logstream=args.logfile):
        kevlar.print_augmented_fastx(contig, outstream)
