#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import khmer
from khmer import khmer_args


def subparser(subparsers):
    """Define the `kevlar filter` command-line interface."""

    desc = """\
    Discard k-mers and reads that are contaminant in origin or whose abundances
    were inflated at the `kevlar count` or `kevlar novel` stage.
    """

    subparser = subparsers.add_parser('filter', description=desc)
    subparser.add_argument(
        '-x', '--ctrl-max', metavar='X', type=int, default=1,
        help='k-mers with abund > X in any control sample are uninteresting; '
        'default is X=1'
    )
    subparser.add_argument(
        '-y', '--case-min', metavar='Y', type=int, default=6,
        help='k-mers with abund < Y in any case sample are uninteresting; '
        'default is Y=6'
    )
    subparser.add_argument(
        '-o', '--out', metavar='FILE', help='output file; default is terminal '
        '(stdout)'
    )
    subparser.add_argument(
        'counts', help='table of k-mer counts re-computed solely from '
        'putatively novel reads'
    )
    subparser.add_argument(
        'augfastq', help='putatively novel reads in augmented Fastq format'
    )
