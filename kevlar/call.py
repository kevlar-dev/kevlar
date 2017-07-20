#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import sys
import khmer
import kevlar


class KevlarAmbiguousTargetError(ValueError):
    pass


def main(args):
    outstream = kevlar.open(args.out, 'w')

    queryseqs = [record for record in khmer.ReadParser(args.queryseq)]
    targetseqs = [record for record in khmer.ReadParser(args.targetseq)]
    nquerys = len(queryseqs)
    ntargets = len(targetseqs)
    if ntargets != 1:
        message = 'found {:d} target sequences, expected 1'.format(ntargets)
        raise KevlarAmbiguousTargetError(message)

    queryseqs = sorted(queryseqs, reverse=True, key=len)

    target = targetseqs[0].sequence
    query = queryseqs[0].sequence
    cigar = kevlar.align(target, query, args.match, args.mismatch, args.open,
                         args.extend)
    print(cigar, file=outstream)
