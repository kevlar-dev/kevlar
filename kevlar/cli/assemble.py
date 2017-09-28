#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import khmer


def subparser(subparsers):
    """Define the `kevlar assemble` command-line interface."""

    desc = "Use a simple greedy algorithm to assemble a single variant's reads"

    subparser = subparsers.add_parser('assemble', description=desc,
                                      add_help=False)

    default_args = subparser.add_argument_group(
        'Default mode',
        'By default, `kevlar assemble` uses a home-grown implementation of a '
        'greedy assembly algorithm. Reads that share interesting k-mers are '
        'merged one by one, starting with the pair of reads with the largest '
        'overlap.'
    )
    default_args.add_argument('-d', '--debug', action='store_true',
                              help='print debugging output')
    default_args.add_argument('--gml', metavar='FILE',
                              help='write graph to .gml file')
    default_args.add_argument('-n', '--min-abund', type=int, metavar='N',
                              default=2, help='discard interesting k-mers '
                              'that occur fewer than N times')
    default_args.add_argument('-x', '--max-abund', type=int, metavar='X',
                              default=500, help='discard interesting k-mers '
                              'that occur more than X times')

    jca_args = subparser.add_argument_group(
        'Junction count assembly mode',
        'Alternatively, `kevlar assemble` can uses a junction count assembler '
        'available from the khmer library.'
    )
    jca_args.add_argument('--jca', action='store_true', help="use khmer's "
                          "junction count assembler instead of kevlar's "
                          "greedy assembler")
    jca_args.add_argument('--ignore', metavar='KM', nargs='+',
                          help='ignore the specified k-mer(s)')
    jca_args.add_argument('--collapse', action='store_true', help='collapse '
                          'contigs that are contained within other contigs')
    jca_args.add_argument('-M', '--memory', default='1e6', metavar='MEM',
                          type=khmer.khmer_args.memory_setting,
                          help='memory to allocate for assembly graph; '
                          'default is 1M')
    jca_args.add_argument('--max-fpr', type=float, metavar='FPR',
                          default=0.001, help='terminate if the expected '
                          'false positive rate of the assembly graph is '
                          'higher than the specified FPR; default is 0.001')

    misc_args = subparser.add_argument_group('Miscellaneous settings')
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')

    subparser.add_argument('augfastq', help='annotated reads in augmented '
                           'Fastq format')
