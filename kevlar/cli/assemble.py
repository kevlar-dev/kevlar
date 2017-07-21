#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    """Define the `kevlar assemble` command-line interface."""

    desc = "Use a simple greedy algorithm to assemble a single variant's reads"

    subparser = subparsers.add_parser('assemble', description=desc)
    subparser.add_argument('-d', '--debug', action='store_true',
                           help='print debugging output')
    subparser.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('--gml', metavar='FILE',
                           help='write graph to .gml file')
    subparser.add_argument('-n', '--min-abund', type=int, metavar='N',
                           default=2, help='discard interesting k-mers that '
                           'occur fewer than N times')
    subparser.add_argument('-x', '--max-abund', type=int, metavar='X',
                           default=500, help='discard interesting k-mers that '
                           'occur more than X times')
    subparser.add_argument('augfastq', help='annotated reads in augmented '
                           'Fastq format')
