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

    refr_args = subparser.add_argument_group(
        'Reference genome',
        'A reference genome is not required, but if available can be supplied '
        'to discard "interesting" k-mers that turn out not to be novel.'
    )
    refr_args.add_argument('--refr', metavar='FILE', type=str, default=None,
                           help='reference genome in Fasta/Fastq format')
    refr_args.add_argument('--refr-memory', metavar='MEM', default='1e6',
                           type=khmer_args.memory_setting,
                           help='memory to allocate for storing the reference '
                           'genome; default is 1M')
    refr_args.add_argument('--refr-max-fpr', type=float, metavar='FPR',
                           default=0.001, help='terminate if the expected '
                           'false positive rate is higher than the specified '
                           'FPR; default is 0.001')

    contam_args = subparser.add_argument_group(
        'Screening for known contaminants',
        'It may not be possible to anticipate every possible contaminant that '
        'may be present in a sample, but if there is a common or known set of '
        'contaminants these can be filtered out at this step. Unknown '
        'contaminants will have to be filtered out after read partitioning.'
    )
    contam_args.add_argument('--contam', metavar='FILE', type=str,
                             default=None, help='database of contaminant '
                             'sequences in Fasta/Fastq format')
    contam_args.add_argument('--contam-memory', metavar='MEM', default='1e6',
                             type=khmer_args.memory_setting,
                             help='memory to allocate for storing contaminant '
                             'sequences; default is 1M')
    contam_args.add_argument('--contam-max-fpr', type=float, metavar='FPR',
                             default=0.001, help='terminate if the expected '
                             'false positive rate is higher than the specified'
                             ' FPR; default is 0.001')

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
    filter_args.add_argument('--abund-memory', metavar='MEM', default='1e6',
                             type=khmer_args.memory_setting,
                             help='memory to allocate for re-calculating '
                             'abundance of interesting k-mers; default is 1M')
    filter_args.add_argument('--abund-max-fpr', type=float, metavar='FPR',
                             default=0.001, help='terminate if the expected '
                             'false positive rate is higher than the specified'
                             ' FPR; default is 0.001')
    filter_args.add_argument('--min-abund', type=int, default=5, metavar='Y',
                             help='minimum abundance required to call a '
                             'k-mer novel; should be the same value used for '
                             '--case_min in `kevlar novel`; default is 5')
    filter_args.add_argument('--skip2', default=False, action='store_true',
                             help='skip the second pass over the reads that '
                             'recalculates abundance after reference and '
                             'contaminant k-mers are discarded')
    filter_args.add_argument('--ignore', metavar='KM', nargs='+',
                             help='ignore the specified k-mer(s)')

    misc_args = subparser.add_argument_group(
        'Miscellaneous settings'
    )
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    misc_args.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    misc_args.add_argument('--aug-out', metavar='FILE',
                           help='optional augmented Fastq output')

    subparser.add_argument('augfastq', nargs='+', help='one or more files in '
                           '"augmented" Fastq format (a la `kevlar novel` '
                           'output)')
