#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import argparse
import sys
import kevlar

def subparser(parser):
    subparsers = parser.add_subparsers(dest='cmd')
    subparser = subparsers.add_parser('dump')
    subparser.add_argument('--infmt', metavar='FMT', choices=['bam, fastq'],
                           help='format of input file; "bam" or "fastq"')
    subparser.add_argument('--out', metavar='FILE',
                           type=argparse.FileType('w'),
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('refr', help='reference sequence in Fasta format')
    subparser.add_argument('reads', help='reads in Fastq or BAM format')


def main(args, log=sys.stderr):
    print('[kevlar::dump] Loading reference sequence', file=log)
    seqs = kevlar.fasta.parse_seq_dict(args.refr)
