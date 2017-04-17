#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict
import argparse
import itertools
import sys
import networkx
import khmer
import kevlar


def subparser(subparsers):
    subparser = subparsers.add_parser('assemble')
    subparser.add_argument('-o', '--out', metavar='FILE', default=sys.stdout,
                           type=argparse.FileType('w'),
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('augfastq', type=argparse.FileType('r'),
                           help='annotated reads in augmented Fastq format')


def calc_offset(read1, read2, minkmer):
    maxkmer = kevlar.revcommin(minkmer)
    kmer1 = [k for k in read1.ikmers
             if kevlar.same_seq(k.sequence, minkmer, maxkmer)]
    kmer2 = [k for k in read2.ikmers
             if kevlar.same_seq(k.sequence, minkmer, maxkmer)]
    ksize = len(kmer1.sequence)

    pos1 = kmer1.offset
    pos2 = kmer2.offset
    sameorient = True
    if kmer1.sequence != kmer2.sequence:
        assert kmer1.sequence == kevlar.revcom(kmer2.sequence)
        sameorient = False
        pos2 = len(read2.sequence) - (kmer2.offset + ksize)

    tail, head = read1, read2
    offset = pos1 - pos2
    if pos2 > pos1:
        tail, head = read2, read1
        offset = pos2 - pos1

    return tail, head, offset


def main(args):
    reads = defaultdict(set)  # key: k-mer (min repr), value: set of reads
    for n, record in enumerate(kevlar.parse_augmented_fastq(args.augfastq), 1):
        if n % 10000 == 0:
            print('[kevlar::assemble]    loaded {:d} reads'.format(n),
                  file=args.logfile)
        for kmer in record.ikmers:
            kmerseq = kevlar.revcommin(kmer.sequence)
            reads[kmerseq].add(record)

    graph = networkx.Graph()  # DiGraph?
    for minkmer in reads:
        readset = reads[minkmer]
        assert len(readset) > 1
        for read1, read2 in itertools.combinations(readset, 2):
            tail, head, offset = calc_offset(read1, read2, minkmer)
            if tail in graph and head in graph[tail]:
                assert graph[tail][head]['offset'] == offset
            else:
                graph.add_edge(tail, head, offset=offset)

    for n, cc in enumerate(networkx.connected_components(graph)):
        print('CC', n, len(cc), sep='\t')
