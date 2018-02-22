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


def augment(augseqstream, nakedseqstream):
    """
    Augment an unannotated stream of sequences.

    - `augseqstream`: a stream of sequences annotated with k-mers of interest
    - `nakedseqstream`: a stream of unannotated sequences, to be augmented with
      k-mers of interest from `augseqstream`
    """
    ksize = None
    ikmers = dict()
    mateseqs = set()
    for record in augseqstream:
        for ikmer in record.ikmers:
            ikmers[ikmer.sequence] = ikmer.abund
            ikmers[kevlar.revcom(ikmer.sequence)] = ikmer.abund
            ksize = len(ikmer.sequence)
        if hasattr(record, 'mateseqs'):
            mateseqs.update(record.mateseqs)
    mateseqs = sorted(mateseqs)

    for record in nakedseqstream:
        newikmers = list()
        numkmers = len(record.sequence) - ksize + 1
        for offset in range(numkmers):
            kmer = record.sequence[offset:offset+ksize]
            if kmer in ikmers:
                ikmer = kevlar.KmerOfInterest(kmer, offset, ikmers[kmer])
                newikmers.append(ikmer)
        if hasattr(record, 'quality'):
            newrecord = screed.Record(
                name=record.name, sequence=record.sequence, ikmers=newikmers,
                quality=record.quality, mateseqs=mateseqs
            )
        else:
            newrecord = screed.Record(
                name=record.name, sequence=record.sequence, ikmers=newikmers,
                mateseqs=mateseqs
            )
        yield newrecord


def main(args):
    augfh = kevlar.open(args.augseqs, 'r')
    augseqs = kevlar.parse_augmented_fastx(augfh)
    nakedseqs = screed.open(args.seqs)
    outstream = kevlar.open(args.out, 'w')
    for record in augment(augseqs, nakedseqs):
        kevlar.print_augmented_fastx(record, outstream)
