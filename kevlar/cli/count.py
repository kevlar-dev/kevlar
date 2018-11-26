#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
import textwrap
import khmer
from khmer import khmer_args


def subparser(subparsers):
    """Define the `kevlar count` command-line interface."""

    desc = """\
    Compute k-mer abundances for the provided sample. Supports k-mer banding:
    see http://kevlar.readthedocs.io/en/latest/banding.html for more details.
    """
    desc = textwrap.dedent(desc)

    epilog = """\
    Example::

        kevlar count --memory 500M case1.ct case1-reads.fastq

    Example::

        kevlar count --ksize 25 --memory 12G --max-fpr 0.01 --threads 8 \\
            proband.counttable \\
            proband-R1.fq.gz proband-R2.fq.gz proband-unpaired.fq.gz"""
    epilog = textwrap.dedent(epilog)

    subparser = subparsers.add_parser(
        'count', description=desc, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparser.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    subparser.add_argument('-c', '--counter-size', type=int, choices=(1, 4, 8),
                           metavar='C', default=8, help='number of bits to '
                           'allocate for counting each k-mer; options are 1 '
                           '(max count: 1), 4 (max count: 15), and 8 (max '
                           'count: 255); default is 8')
    subparser.add_argument('-M', '--memory', type=khmer_args.memory_setting,
                           default=1e6, metavar='MEM',
                           help='memory to allocate for the count table')
    subparser.add_argument('--max-fpr', type=float, default=0.2,
                           metavar='FPR', help='terminate if the estimated '
                           'false positive rate for any sample is higher than '
                           '"FPR"; default is 0.2')
    subparser.add_argument('--mask', metavar='MSK', help='counttable or '
                           'nodetable of k-mers to ignore when counting '
                           'k-mers')
    subparser.add_argument('--count-masked', action='store_true',
                           help='by default, when a mask is provided k-mers '
                           'in the mask are ignored; this setting inverts the '
                           'behavior so that only k-mers in the mask are '
                           'counted')
    subparser.add_argument('--num-bands', type=int, metavar='N', default=None,
                           help='number of bands into which to divide the '
                           'hashed k-mer space')
    subparser.add_argument('--band', type=int, metavar='I', default=None,
                           help='a number between 1 and N (inclusive) '
                           'indicating the band to be processed')
    subparser.add_argument('-t', '--threads', type=int, default=1, metavar='T',
                           help='number of threads to use for file processing;'
                           ' default is 1')
    subparser.add_argument('counttable', type=str, help='name of the file to '
                           'which the output (a k-mer count table) will be '
                           'written; the suffix ".counttable" will be applied '
                           'if the provided file name does not end in ".ct" '
                           'or ".counttable"')
    subparser.add_argument('seqfile', type=str, nargs='+', help='input files '
                           'in Fastq/Fasta format')
