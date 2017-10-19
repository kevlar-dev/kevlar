#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
import sys
import screed
import kevlar


def reaugment(augreadstream, nakedreadstream):
    augreads = dict()
    for record in augreadstream:
        augreads[record.name] = record

    for record in nakedreadstream:
        augrecord = augreads[record.name]
        record.ikmers = augrecord.ikmers
        if len(record.sequence) != len(augrecord.sequence):
            assert len(record.sequence) < len(augrecord.sequence)
            ikmers = list()
            for kmer in augrecord.ikmers:
                stillthere = (
                    kmer.sequence in record.sequence or
                    kevlar.revcom(kmer.sequence) in record.sequence
                )
                if stillthere:
                    ikmers.append(kmer)
            if len(ikmers) == 0:
                continue
            record.ikmers = ikmers
        yield record


def main(args):
    augfh = kevlar.open(args.augfastq, 'r')
    augreads = kevlar.parse_augmented_fastx(augfh)
    nakedreads = screed.open(args.fastq)
    outstream = kevlar.open(args.out, 'w')
    for record in reaugment(augreads, nakedreads):
        kevlar.print_augmented_fastx(record, outstream)
