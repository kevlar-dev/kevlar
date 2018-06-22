#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import khmer


def subparser(subparsers):
    """Define the `kevlar assemble` command-line interface."""

    desc = 'Assemble reads into contigs representing putative variants'
    subparser = subparsers.add_parser('assemble', description=desc)
    subparser.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('augfastq', help='annotated reads in augmented '
                           'Fastq/Fasta format')
