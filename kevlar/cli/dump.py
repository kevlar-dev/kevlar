#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import khmer
from khmer import khmer_args


def subparser(subparsers):
    """Define the `kevlar dump` command-line interface."""

    desc = 'Discard reads that align perfectly to the reference genome.'

    subparser = subparsers.add_parser('dump', description=desc)
    subparser.add_argument('--seqid', metavar='SEQ',
                           help='dump reads not mapped to SEQ')
    subparser.add_argument('--genomemask', metavar='FILE', help='dump reads '
                           'with median k-mer abundance >= 1 in the specified '
                           'genome; if both --seqid and --genomemask are '
                           'declared, reads passing either filter will be '
                           'kept')
    subparser.add_argument('--maskmemory', metavar='SIZE', default=2e9,
                           type=khmer_args.memory_setting,
                           help='memory to be occupied by genome mask; default'
                           ' is 2G')
    subparser.add_argument('--mask-k', metavar='K', default=31, type=int,
                           help='k size for genome mask')
    subparser.add_argument('--out', metavar='FILE', help='output file; default'
                           ' is terminal (stdout)')
    subparser.add_argument('refr', help='reference sequence in Fasta format')
    subparser.add_argument('reads', help='read alignments in BAM format')
