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
    subparser.add_argument('-o', '--out', metavar='FILE', default=sys.stdout,
                           type=argparse.FileType('w'),
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('augfastq', type=argparse.FileType('r'),
                           help='original augmented Fastq file')
    subparser.add_argument('fastq', help='processed Fastq file to re-annotate')


def main(args):
    reads = dict()
    for record in kevlar.parse_augmented_fastq(args.augfastq):
        reads[record.name] = record

    reader = khmer.ReadParser(args.fastq)
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
        kevlar.print_augmented_fastq(augrecord, args.out)
