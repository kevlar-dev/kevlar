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
import sys
import khmer
import kevlar


def subparser(subparsers):
    subparser = subparsers.add_parser('reaugment')
    subparser.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('augfastq', help='original augmented Fastq file')
    subparser.add_argument('fastq', help='processed Fastq file to re-annotate')


def main(args):
    reads = dict()
    instream = kevlar.open(args.augfastq, 'r')
    for record in kevlar.parse_augmented_fastq(instream):
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
        kevlar.print_augmented_fastq(augrecord, outstream)
