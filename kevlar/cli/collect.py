#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
import khmer
from khmer import khmer_args


def subparser(subparsers):
    """Define the `kevlar collect` command-line interface."""

    subparser = subparsers.add_parser('collect')
    subparser.add_argument('-d', '--debug', action='store_true',
                           help='print debugging output')
    subparser.add_argument('-M', '--memory', default='1e6', metavar='MEM',
                           type=khmer_args.memory_setting,
                           help='memory to allocate for recalculating '
                           'abundances of novel k-mers; default is 1M')
    subparser.add_argument('--minabund', type=int, default=5, metavar='Y',
                           help='minimum case abundance required to call a '
                           'k-mer novel; used to filter out k-mers with '
                           'inflated abundances; should equal the value of '
                           '--case_min used in "kevlar novel"; default is 5')
    subparser.add_argument('--max-fpr', type=float, metavar='FPR',
                           default=0.001, help='terminate if the expected '
                           'false positive rate is higher than the specified '
                           'FPR; default is 0.001')
    subparser.add_argument('--refr', metavar='FILE', type=str, default=None,
                           help='reference genome in Fasta/Fastq format; any '
                           'k-mers designated as "interesting" by "kevlar '
                           'novel" are ignored if they are present in the '
                           'reference genome')
    subparser.add_argument('--refr-memory', metavar='MEM', default='1e6',
                           type=khmer_args.memory_setting,
                           help='memory to allocate for storing the reference '
                           'genome; default is 1M')
    subparser.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    subparser.add_argument('--ignore', metavar='KMER', nargs='+',
                           help='ignore the specified k-mer(s)')
    subparser.add_argument('-o', '--out', type=argparse.FileType('w'),
                           metavar='OUT',
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('--collapse', action='store_true', help='collapse '
                           'linear paths contained in other linear paths')
    subparser.add_argument('novel_output', nargs='+', help='one or more output'
                           ' files from the "kevlar novel" command')
