#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import khmer
import kevlar


def call(targetlist, querylist, match=1, mismatch=2, gapopen=5, gapextend=0):
    """
    Wrap the `kevlar call` procedure as a generator function.

    Input is the following.
    - an iterable containing one or more target sequences from the reference
      genome, stored as khmer or screed sequence records
    - an iterable containing one or more contigs assembled by kevlar, stored as
      khmer or screed sequence records
    - alignment match score (integer)
    - alignment mismatch penalty (integer)
    - alignment gap open penalty (integer)
    - alignment gap extension penalty (integer)

    The function yields tuples of target sequence name, query sequence name,
    and alignment CIGAR string
    """
    for target in sorted(targetlist, key=lambda record: record.name):
        for query in sorted(querylist, reverse=True, key=len):
            cigar = kevlar.align(target.sequence, query.sequence, match,
                                 mismatch, gapopen, gapextend)
            yield target.name, query.name, cigar


def main(args):
    outstream = kevlar.open(args.out, 'w')
    qinstream = kevlar.parse_augmented_fastx(kevlar.open(args.queryseq, 'r'))
    queryseqs = [record for record in qinstream]
    targetseqs = [record for record in khmer.ReadParser(args.targetseq)]
    caller = call(
        targetseqs, queryseqs,
        args.match, args.mismatch, args.open, args.extend
    )
    for target, query, cigar in caller:
        print(target, query, cigar, sep='\t', file=outstream)
