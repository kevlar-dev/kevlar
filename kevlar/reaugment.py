#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict
import sys
import khmer
import kevlar


def main(args):
    reads = dict()
    instream = kevlar.open(args.augfastq, 'r')
    for record in kevlar.parse_augmented_fastx(instream):
        reads[record.name] = record

    reader = khmer.ReadParser(args.fastq)
    outstream = kevlar.open(args.out, 'w')
    for read in reader:
        augrecord = reads[read.name]
        if len(read.sequence) < len(augrecord.sequence):
            ikmers = list()
            for kmer in augrecord.ikmers:
                stillthere = (
                    kmer.sequence in read.sequence or
                    kevlar.revcom(kmer.sequence) in read.sequence
                )
                if stillthere:
                    ikmers.append(kmer)
            if len(ikmers) == 0:
                continue
            augrecord.ikmers = ikmers
        kevlar.print_augmented_fastx(augrecord, outstream)
