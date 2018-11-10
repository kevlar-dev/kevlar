#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
import sys


def augment(augseqstream, nakedseqstream, collapsemates=False, upint=10000):
    """
    Augment an unannotated stream of sequences.

    - `augseqstream`: a stream of sequences annotated with k-mers of interest
    - `nakedseqstream`: a stream of unannotated sequences, to be augmented with
      k-mers of interest from `augseqstream`
    """
    ksize = None
    ikmers = dict()
    mateseqs = dict()
    for n, record in enumerate(augseqstream):
        if n > 0 and n % upint == 0:
            print('[kevlar::augment] processed', n, 'input reads',
                  file=sys.stderr)
        for ikmer in record.annotations:
            seq = record.ikmerseq(ikmer)
            ikmers[seq] = ikmer.abund
            ikmers[kevlar.revcom(seq)] = ikmer.abund
            ksize = ikmer.ksize
        assert len(record.mates) in (0, 1)
        if len(record.mates) == 1:
            mateseqs[record.name] = record.mates[0]

    for record in nakedseqstream:
        qual = None
        if hasattr(record, 'quality') and record.quality is not None:
            qual = record.quality
        mates = list()
        if collapsemates:
            mates = sorted(mateseqs.values())
        else:
            if record.name in mateseqs:
                mates.append(mateseqs[record.name])
        newrecord = kevlar.sequence.Record(
            name=record.name, sequence=record.sequence, quality=qual,
            mates=mates,
        )
        numkmers = len(record.sequence) - ksize + 1
        for offset in range(numkmers):
            kmer = record.sequence[offset:offset+ksize]
            if kmer in ikmers:
                abund = ikmers[kmer]
                newrecord.annotate(kmer, offset, abund)
        yield newrecord


def main(args):
    augseqs = kevlar.parse_augmented_fastx(kevlar.open(args.augseqs, 'r'))
    nakedseqs = kevlar.parse_augmented_fastx(kevlar.open(args.seqs, 'r'))
    outstream = kevlar.open(args.out, 'w')
    docollapse = args.collapse_mates
    for record in augment(augseqs, nakedseqs, collapsemates=docollapse):
        kevlar.print_augmented_fastx(record, outstream)
