#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    """Define the `kevlar novel` command-line interface."""

    desc = """\
    Construct a graph to group reads by shared interesting k-mers. Nodes in the
    graph represent reads, and edges between a pair of nodes indicate that the
    two corresponding reads have one or more interesting k-mers in common.
    Connected components in the undirected graph correspond to distinct
    variants (or variant-related breakpoints).
    """

    subparser = subparsers.add_parser('partition', description=desc)
    subparser.add_argument('-d', '--debug', action='store_true',
                           help='print debugging output')
    subparser.add_argument('-s', '--strict', action='store_true',
                           help='require perfect identity between overlapping '
                           'reads for inclusion in the same partition; by '
                           'default, only a shared interesting k-mer is '
                           'required')
    subparser.add_argument('--min-abund', metavar='X', type=int, default=2,
                           help='ignore k-mers with abundance lower than X; '
                           'default is 2')
    subparser.add_argument('--max-abund', metavar='Y', type=int, default=200,
                           help='ignore k-mers with abundance higher than Y; '
                           'default is 200')
    subparser.add_argument('--gml', metavar='FILE',
                           help='write read graph to .gml file')
    subparser.add_argument('outprefix', help='prefix for output files')
    subparser.add_argument('augfastq', help='original augmented Fastq file')
