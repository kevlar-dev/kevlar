#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    """Define the `kevlar split` command-line interface."""

    desc = "Split partitioned reads into N output files."
    subparser = subparsers.add_parser('split', description=desc)
    subparser.add_argument('infile', help='input file; partitioned reads in '
                           'augmented Fastq/Fasta format')
    subparser.add_argument('numfiles', type=int, help='number of output files '
                           'to create')
    subparser.add_argument('base', help='prefix of all output files')
