#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
import sys
import screed
import kevlar


def augment(augreadstream, nakedseqtream):
    ksize = None
    ikmers = dict()
    for record in augreadstream:
        for ikmer in record.ikmers:
            ikmers[ikmer.sequence] = ikmer.abund
            ikmers[kevlar.revcom(ikmer.sequence)] = ikmer.abund
            ksize = len(ikmer.sequence)

    for record in nakedseqtream:
        newikmers = list()
        numkmers = len(record.sequence) - ksize + 1
        for offset in range(numkmers):
            kmer = record.sequence[offset:offset+ksize]
            if kmer in ikmers:
                ikmer = kevlar.KmerOfInterest(kmer, offset, ikmers[kmer])
                newikmers.append(ikmer)
        newrecord = screed.Record(name=record.name, sequence=record.sequence,
                                  ikmers=newikmers)
        yield newrecord


def main(args):
    augfh = kevlar.open(args.augfastq, 'r')
    augreads = kevlar.parse_augmented_fastx(augfh)
    nakedreads = screed.open(args.sequences)
    outstream = kevlar.open(args.out, 'w')
    for record in augment(augreads, nakedreads):
        kevlar.print_augmented_fastx(record, outstream)
