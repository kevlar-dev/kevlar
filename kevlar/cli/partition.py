#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    subparser = subparsers.add_parser('partition')
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
