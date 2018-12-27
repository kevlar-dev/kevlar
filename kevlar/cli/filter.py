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
    were inflated during the preliminary k-mer counting stage.
    """

    subparser = subparsers.add_parser('recount', description=desc)
    subparser.add_argument(
        '-M', '--memory', type=khmer_args.memory_setting, default=1e6,
        metavar='MEM', help='memory to allocate for the k-mer re-counting'
    )
    subparser.add_argument(
        '--max-fpr', type=float, default=0.01, metavar='FPR',
        help='terminate early if the estimated false positive rate for re-'
        'computed k-mer abundances is higher than "FPR"; default is 0.01'
    )
    subparser.add_argument(
        '--mask', metavar='MSK', help='counttable or nodetable of k-mers to '
        'ignore when re-counting k-mers'
    )
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
        'augfastq', help='putatively novel reads in augmented Fastq format'
    )
