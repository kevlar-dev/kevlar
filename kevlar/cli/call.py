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
    """Define the `kevlar call` command-line interface."""

    desc = """\
    Align variant-related reads to the reference genome and call the variant
    from the alignment.
    """

    subparser = subparsers.add_parser('call', description=desc, add_help=False)
    subparser._positionals.title = 'Required inputs'

    score_args = subparser.add_argument_group('Alignment scoring')
    score_args.add_argument('-A', '--match', type=int, default=1, metavar='A',
                            help='match score; default is 1')
    score_args.add_argument('-B', '--mismatch', type=int, default=2,
                            metavar='B', help='mismatch penalty; default is 2')
    score_args.add_argument('-O', '--open', type=int, default=5, metavar='O',
                            help='gap open penalty; default is 5')
    score_args.add_argument('-E', '--extend', type=int, default=0, metavar='E',
                            help='gap extension penalty; default is 0')

    mask_args = subparser.add_argument_group('Mask generation settings')
    mask_args.add_argument('--gen-mask', metavar='FILE', help='generate '
                           'a nodetable containing all k-mers that span any '
                           'variant call')
    mask_args.add_argument('--mask-mem', type=khmer_args.memory_setting,
                           default=1e6, metavar='MEM',
                           help='memory to allocate for the node table')
    mask_args.add_argument('--mask-max-fpr', type=float, default=0.01,
                           metavar='FPR', help='terminate if the estimated '
                           'false positive rate is higher than "FPR"; default '
                           'is 0.01')

    misc_args = subparser.add_argument_group('Miscellaneous settings')
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('-d', '--debug', action='store_true',
                           help='show debugging output')
    misc_args.add_argument('--refr', metavar='FILE',
                           help='reference genome indexed for BWA search; if '
                           'provided, mates of interesting reads will be used '
                           'to diambiguate multi-mapping contigs')
    misc_args.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    misc_args.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    subparser.add_argument('queryseq', help='contigs assembled by "kevlar '
                           'assemble"')
    subparser.add_argument('targetseq', help='region of reference genome '
                           'identified by "kevlar localize"')
