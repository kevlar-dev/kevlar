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
    subparser.add_argument('--out', metavar='FILE', help='output file; default'
                           ' is terminal (stdout)')
    subparser.add_argument('refr', help='reference sequence in Fasta format')
    subparser.add_argument('reads', help='read alignments in BAM format')
