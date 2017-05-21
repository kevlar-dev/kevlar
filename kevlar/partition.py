#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import argparse
import sys
import networkx
import khmer
import kevlar
from kevlar.seqio import load_reads_and_kmers


def subparser(subparsers):
    subparser = subparsers.add_parser('partition')
    subparser.add_argument('-d', '--debug', action='store_true',
                           help='print debugging output')
    subparser.add_argument('-s', '--strict', action='store_true',
                           help='require perfect identity between overlapping '
                           'reads for inclusion in the same partition; by '
                           'default, only a shared interesting k-mer is '
                           'required')
    subparser.add_argument('--gml', metavar='FILE',
                           help='write read graph to .gml file')
    subparser.add_argument('outprefix', help='prefix for output files')
    subparser.add_argument('augfastq', help='original augmented Fastq file')


def main(args):
    debugout = None
    if args.debug:
        debugout = args.logfile

    inputfile = kevlar.open(args.augfastq, 'r')
    reads, kmers = load_reads_and_kmers(inputfile, debugout)
    inputreads = list(reads)
    graph = kevlar.overlap.graph_init(reads, kmers, args.min_abund,
                                      args.max_abund, debugout)
    if args.gml:
        kevlar.to_gml(graph, args.gml, args.logfile)
    kevlar.overlap.write_partitions(graph, reads, args.ccprefix, args.logfile)
