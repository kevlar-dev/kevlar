#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    """Define the `kevlar augment` command-line interface."""

    desc = """\
    Annotate the given sequences with interested k-mers annotated in the given
    augmented Fastq file.
    """

    subparser = subparsers.add_parser('augment', description=desc)
    subparser.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('augfastq', help='original augmented Fastq file')
    subparser.add_argument('sequences', help='sequences to annotate')
