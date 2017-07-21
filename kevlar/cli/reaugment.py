#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    """Define the `kevlar reaugment` command-line interface."""

    desc = """\
    Internally, kevlar uses an "augmented" Fastq format for storing "interesing
    k-mers" and the reads that contain them. Integrating kevlar with other
    tools, however, requires discarding this extra information. The "kevlar
    reaugment" command is used to re-annotate interesting k-mers in Fastq
    files where the information has been stripped.
    """

    subparser = subparsers.add_parser('reaugment', description=desc)
    subparser.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('augfastq', help='original augmented Fastq file')
    subparser.add_argument('fastq', help='processed Fastq file to re-annotate')
