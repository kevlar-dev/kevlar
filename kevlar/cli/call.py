#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


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

    like_args = subparser.add_argument_group('Likelihood filtering')
    like_args.add_argument('--like-filter', action='store_true',
                           help='compute likelihoods for all variant calls')
    like_args.add_argument('--mu', type=float, default=[30.0], nargs='+',
                           metavar='μ', help='mean k-mer abundance (default: '
                           '30.0); can provide a single value for all samples '
                           'or 1 value per sample')
    like_args.add_argument('--sigma', type=float, default=[8.0], nargs='+',
                           metavar='σ', help='k-mer abundance standard '
                           'deviation (default: 8.0); can provide a single '
                           'value for all samples or 1 value per sample')
    like_args.add_argument('--epsilon', type=float, default=[0.01], nargs='+',
                           metavar='ε', help='error rate (default: 0.01); can '
                           'provide a single value for all samples or 1 value '
                           'per sample')
    like_args.add_argument('-r', '--refr', metavar='NT',
                           help='nodetable of reference genome')
    like_args.add_argument('--case', metavar='CT',
                           help='counttable of proband k-mer counts')
    like_args.add_argument('--ctrl', metavar='CT', nargs='+',
                           help='counttables of parent k-mer counts, 1 per '
                           'sample')

    misc_args = subparser.add_argument_group('Miscellaneous settings')
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('--refridx', metavar='FILE',
                           help='reference genome indexed for BWA search; if '
                           'provided, mates of interesting reads will be used '
                           'to diambiguate multi-mapping contigs')
    misc_args.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    misc_args.add_argument('--case-label', metavar='LBL', help='name or label '
                           'for proband/case sample')
    misc_args.add_argument('--ctrl-labels', metavar='LBL', help='comma-'
                           'separated list of names or labels for parent/'
                           'control samples')
    misc_args.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    subparser.add_argument('queryseq', help='contigs assembled by "kevlar '
                           'assemble"')
    subparser.add_argument('targetseq', help='region of reference genome '
                           'identified by "kevlar localize"')
