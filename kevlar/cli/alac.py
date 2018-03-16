#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


def subparser(subparsers):

    desc = """\
    Assemble partitioned reads, localize to the reference genome, align the
    assembled contig to the reference, and call variants.
    """

    subparser = subparsers.add_parser('alac', description=desc, add_help=False)
    subparser._positionals.title = 'Required inputs'

    asmbl_args = subparser.add_argument_group('Read assembly')
    asmbl_args.add_argument('-p', '--part-id', type=str, metavar='ID',
                            help='only process partition "PID" in the input')
    asmbl_args.add_argument('--bigpart', type=int, metavar='N', default=10000,
                            help='do not attempt to assembly any partitions '
                            'with more than N reads (default: 10000)')
    asmbl_args.add_argument('--greedy', action='store_true', help='Use a home-'
                            'grown greedy assembly algorithm instead of the '
                            'default fermi-lite algorithm')
    asmbl_args.add_argument('--fallback', action='store_true', dest='fallback',
                            help='attempt to assemble reads with a home-grown '
                            'greedy assembly algorithm if fermi-lite fails to '
                            'assemble a contig for a partition')

    local_args = subparser.add_argument_group('Target extraction')
    local_args.add_argument('-z', '--seed-size', type=int, default=51,
                            metavar='Z', help='seed size; default is 51')
    local_args.add_argument('-d', '--delta', type=int, default=50, metavar='D',
                            help='retrieve the genomic interval from the '
                            'reference by extending beyond the span of all '
                            'k-mer starting positions by D bp')
    local_args.add_argument('-x', '--max-diff', type=int, metavar='X',
                            default=10000, help='span of all k-mer starting '
                            'positions should not exceed X bp; default is '
                            '10000 (10 kb)')

    score_args = subparser.add_argument_group('Alignment scoring')
    score_args.add_argument('-A', '--match', type=int, default=1, metavar='A',
                            help='match score; default is 1')
    score_args.add_argument('-B', '--mismatch', type=int, default=2,
                            metavar='B', help='mismatch penalty; default is 2')
    score_args.add_argument('-O', '--open', type=int, default=5, metavar='O',
                            help='gap open penalty; default is 5')
    score_args.add_argument('-E', '--extend', type=int, default=0, metavar='E',
                            help='gap extension penalty; default is 0')

    misc_args = subparser.add_argument_group('Miscellaneous settings')
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    misc_args.add_argument('-i', '--min-ikmers', metavar='I', type=int,
                           default=None, help='do not report calls that a '
                           'supported by fewer than `I` interesting k-mers')
    misc_args.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    subparser.add_argument('infile', help='partitioned reads in augmented '
                           'Fastq format')
    subparser.add_argument('refr', help='reference genome in Fasta format '
                           '(indexed for bwa search)')
