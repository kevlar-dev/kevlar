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
    were inflated at the "kevlar count" or "kevlar novel" stage.
    """

    subparser = subparsers.add_parser('filter', description=desc,
                                      add_help=False)
    subparser._positionals.title = 'Required inputs'

    mask_args = subparser.add_argument_group(
        'Reference and contaminant masking',
        'Two classes of "interesting k-mers" (ikmers) need to be filtered out '
        'when identifying novel germline mutations: those in the reference '
        'genome, and those of contaminant origin. When available, provide '
        'reference genome assemblies and contaminant databases to "mask" out '
        'these k-mers and remove their "interesting" classification.'
    )
    mask_args.add_argument('--mask', metavar='FA', type=str, default=None,
                           nargs='+', help='sequences to mask (reference '
                           'genomes, contaminants); can provide as one or more'
                           ' Fasta/Fastq files or as a single pre-computed '
                           'nodetable file; see `--save-mask` option')
    mask_args.add_argument('--mask-memory', metavar='MEM', default='1e9',
                           type=khmer_args.memory_setting,
                           help='memory to allocate for storing the mask; '
                           'default is 1G')
    mask_args.add_argument('--mask-max-fpr', type=float, metavar='FPR',
                           default=0.001, help='terminate if the expected '
                           'false positive rate is higher than the specified '
                           'FPR; default is 0.001')
    mask_args.add_argument('--save-mask', type=str, metavar='FILE',
                           default=None, help='save mask nodetable to the '
                           'specified FILE to save time with future runs')

    filter_args = subparser.add_argument_group(
        'Filtering k-mers',
        'Memory constraints often require running `kevlar novel` with false '
        'positive rates (FPRs) in the 0.1 - 0.2 range, resulting in some '
        'k-mers reported with highly inflated abundances. This script handles '
        'a much smaller amount of data, and in limited memory can achieve a '
        'much lower FPR, compute exact k-mer abundances, and discard '
        '"interesting" k-mers whose abundances were incorrectly reported '
        'previously.'
    )
    filter_args.add_argument(
        '--abund-memory', metavar='MEM', default='1e6',
        type=khmer_args.memory_setting, help='memory to allocate for re-'
        'calculating abundance of interesting k-mers; default is 1M'
    )
    filter_args.add_argument(
        '--abund-max-fpr', type=float, metavar='FPR', default=0.001,
        help='terminate if the expected false positive rate is higher than '
        'the specified FPR; default is 0.001'
    )
    filter_args.add_argument(
        '-x', '--ctrl-max', metavar='X', type=int, default=1,
        help='k-mers with abund > X in any control sample are uninteresting; '
        'default is X=1'
    )
    filter_args.add_argument(
        '-y', '--case-min', metavar='Y', type=int, default=6,
        help='k-mers with abund < Y in any case sample are uninteresting; '
        'default is Y=6'
    )
    filter_args.add_argument(
        '--ignore', metavar='KM', nargs='+',
        help='ignore the specified k-mer(s)'
    )

    misc_args = subparser.add_argument_group(
        'Miscellaneous settings'
    )
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    misc_args.add_argument('-o', '--out', metavar='FILE',
                           help='output file (in augmented Fastq format); '
                           'default is terminal (stdout)')

    subparser.add_argument('augfastq', nargs='+', help='one or more files in '
                           '"augmented" Fastq format (a la `kevlar novel` '
                           'output)')
