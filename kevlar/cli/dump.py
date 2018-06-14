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
    subparser.add_argument('-o', '--out', metavar='FILE', help='output file; '
                           'default is terminal (stdout)')
    subparser.add_argument('-r', '--refr', metavar='REFR', help='reference '
                           'sequence in Fasta format')
    subparser.add_argument('--pair-mode', choices=('split', 'keep', 'drop'),
                           default='split', metavar='PM',
                           help='how to handle situations where only one read '
                           'in a pair qualifies to be discarded; in "split" '
                           'mode, discard one read and keep the other; in '
                           '"keep" mode, keep the entire pair; in "drop" '
                           'mode, discard the entire pair; default is "split"')
    subparser.add_argument('reads', help='read alignments in BAM format')
