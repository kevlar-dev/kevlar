#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
import sys


def assemble_fml_asm(readstream, logstream=sys.stderr):
    reads = [r for r in readstream]
    assembler = kevlar.assembly.fml_asm(reads)
    for n, contig in enumerate(assembler, 1):
        name = 'contig{:d}'.format(n)
        record = kevlar.sequence.Record(name=name, sequence=contig)
        yield next(kevlar.augment.augment(reads, [record]))


def main(args):
    reads = kevlar.parse_augmented_fastx(kevlar.open(args.augfastq, 'r'))
    outstream = kevlar.open(args.out, 'w')
    for contig in assemble_fml_asm(reads):
        kevlar.print_augmented_fastx(contig, outstream)
