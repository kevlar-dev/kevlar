#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import re


def subparser(subparsers):
    """Define the `kevlar localize` command-line interface."""

    desc = """\
    For each partition, compute the reference target sequence to use for
    variant calling. NOTE: this command relies on the `bwa` program being in
    the PATH environmental variable.
    """

    subparser = subparsers.add_parser('localize', description=desc)
    subparser.add_argument('-d', '--delta', type=int, metavar='D',
                           default=50, help='retrieve the genomic interval '
                           'from the reference by extending beyond the span '
                           'of all k-mer starting positions by D bp')
    subparser.add_argument('-p', '--part-id', type=str, metavar='ID',
                           help='only localize partition "ID" in the input')
    subparser.add_argument('-o', '--out', metavar='FILE', default='-',
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('-z', '--seed-size', type=int, metavar='Z',
                           default=51, help='seed size; default is 51')
    subparser.add_argument('-x', '--max-diff', type=int, metavar='X',
                           default=None, help='split and report multiple '
                           'reference targets if the distance between two '
                           'seed matches is > X; by default, X is set '
                           'dynamically for each partition and is equal to 3 '
                           'times the length of the longest contig in the '
                           'partition; each resulting bin specifies a '
                           'reference target sequence to which assembled '
                           'contigs will subsequently be aligned')
    subparser.add_argument('--include', metavar='REGEX', type=str,
                           help='discard alignments to any chromosomes whose '
                           'sequence IDs do not match the given pattern')
    subparser.add_argument('--exclude', metavar='REGEX', type=str,
                           help='discard alignments to any chromosomes whose '
                           'sequence IDs match the given pattern')
    subparser.add_argument('contigs', help='assembled reads in Fasta format')
    subparser.add_argument('refr', help='BWA indexed reference genome')
