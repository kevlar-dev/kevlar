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
import re


def subparser(subparsers):

    desc = """\
    Assemble partitioned reads, localize to the reference genome, align the
    assembled contig to the reference, and call variants.
    """

    subparser = subparsers.add_parser('alac', description=desc, add_help=False)
    subparser._positionals.title = 'Required inputs'

    asmbl_args = subparser.add_argument_group('Read assembly')
    asmbl_args.add_argument('-p', '--part-id', type=str, metavar='ID',
                            help='only process partition "ID" in the input')
    asmbl_args.add_argument('--max-reads', type=int, metavar='N',
                            default=10000, help='do not attempt to assemble '
                            'any partitions with more than N reads (default: '
                            '10000)')

    local_args = subparser.add_argument_group('Target extraction')
    local_args.add_argument('-z', '--seed-size', type=int, default=51,
                            metavar='Z', help='seed size; default is 51')
    local_args.add_argument('-d', '--delta', type=int, default=50, metavar='D',
                            help='retrieve the genomic interval from the '
                            'reference by extending beyond the span of all '
                            'k-mer starting positions by D bp')
    local_args.add_argument('-x', '--max-diff', type=int, metavar='X',
                            default=None, help='split and report multiple '
                            'reference targets if the distance between two '
                            'seed matches is > X; by default, X is set '
                            'dynamically for each partition and is equal to 3 '
                            'times the length of the longest contig in the '
                            'partition; each resulting bin specifies a '
                            'reference target sequence to which assembled '
                            'contigs will subsequently be aligned')
    local_args.add_argument('--include', metavar='REGEX', type=str,
                            help='discard alignments to any chromosomes whose '
                            'sequence IDs do not match the given pattern')
    local_args.add_argument('--exclude', metavar='REGEX', type=str,
                            help='discard alignments to any chromosomes whose '
                            'sequence IDs match the given pattern')

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
    misc_args.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    misc_args.add_argument('-i', '--min-ikmers', metavar='I', type=int,
                           default=None, help='do not report calls that a '
                           'supported by fewer than `I` interesting k-mers')
    misc_args.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    misc_args.add_argument('-t', '--threads', type=int, default=1, metavar='T',
                           help='process T partitions at a time using T '
                           'threads')
    subparser.add_argument('infile', help='partitioned reads in augmented '
                           'Fastq format')
    subparser.add_argument('refr', help='reference genome in Fasta format '
                           '(indexed for bwa search)')
