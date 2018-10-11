#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import re


def subparser(subparsers):
    """Define the `kevlar localize` command-line interface."""

    desc = """\
    Given a reference genome and a contig (or set of contigs) assembled from
    variant-related reads, retrieve the portion of the reference genome
    corresponding to the variant. NOTE: this command relies on the `bwa`
    program being in the PATH environmental variable.
    """

    subparser = subparsers.add_parser('localize', description=desc)
    subparser.add_argument('-d', '--delta', type=int, metavar='D',
                           default=50, help='retrieve the genomic interval '
                           'from the reference by extending beyond the span '
                           'of all k-mer starting positions by D bp')
    subparser.add_argument('-o', '--out', metavar='FILE', default='-',
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('-z', '--seed-size', type=int, metavar='Z',
                           default=51, help='seed size; default is 51')
    subparser.add_argument('-x', '--max-diff', type=int, metavar='X',
                           default=None, help='split and report multiple '
                           'reference targets if the distance between two '
                           'seed matches is > X; by default, X is 3 times the '
                           'length of the longest contig')
    subparser.add_argument('--include', metavar='REGEX', type=re.escape,
                           help='discard alignments to any chromosomes whose '
                           'sequence IDs do not match the given pattern')
    subparser.add_argument('--exclude', metavar='REGEX', type=re.escape,
                           help='discard alignments to any chromosomes whose '
                           'sequence IDs match the given pattern')
    subparser.add_argument('contigs', help='assembled reads in Fasta format')
    subparser.add_argument('refr', help='BWA indexed reference genome')

