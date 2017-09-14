#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import re
import sys
import khmer
import kevlar


def make_call(target, query, cigar):
    snvmatch = re.search('(\d+)D(\d+)M(\d+)D', cigar)
    if snvmatch:
        offset = int(snvmatch.group(1))
        length = int(snvmatch.group(2))
        t = target[offset:offset+length]
        assert len(t) == length
        q = query[:length]
        assert len(q) == length
        diffs = [i for i in range(length) if t[i] != q[i]]
        if len(diffs) == 1:
            return offset + diffs[0]


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
            varcall = make_call(target.sequence, query.sequence, cigar)
            yield target.name, query.name, cigar, varcall


def main(args):
    outstream = kevlar.open(args.out, 'w')
    qinstream = kevlar.parse_augmented_fastx(kevlar.open(args.queryseq, 'r'))
    queryseqs = [record for record in qinstream]
    targetseqs = [record for record in khmer.ReadParser(args.targetseq)]
    caller = call(
        targetseqs, queryseqs,
        args.match, args.mismatch, args.open, args.extend
    )
    for target, query, cigar, varcall in caller:
        print(target, query, cigar, varcall, sep='\t', file=outstream)
