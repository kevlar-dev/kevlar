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

    desc = "Compute k-mer abundances for the provided samples"
    epilog = """\
    Example::

        kevlar count --ksize 31 --memory 4G --mem-frac 0.25 \\
            --case proband.counttable proband-reads.fq.gz \\
            --control father.counttable father-r1.fq.gz father-r2.fq.gz \\
            --control mother.counttable mother-reads.fq.gz

    Example::

        kevlar count --ksize 25 --memory 500M \\
            --case case_x.ct case_x.fq \\
            --case case_y.ct case_y.fq \\
            --control control1.ct control1a.fq control1b.fq \\
            --control control2.ct control2a.fq control2b.fq"""
    epilog = textwrap.dedent(epilog)

    subparser = subparsers.add_parser(
        'count', description=desc, epilog=epilog, add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    samp_args = subparser.add_argument_group('Case/control config')
    samp_args.add_argument(
        '--case', metavar='F', nargs='+', required=True, action='append',
        help='one or more FASTA/FASTQ files containing reads from a case '
        'sample; can be declared multiple times corresponding to multiple '
        'case samples; see examples below'
    )
    samp_args.add_argument(
        '--control', metavar='F', nargs='+', action='append',
        help='one or more FASTA/FASTQ files containing reads from a control '
        'sample; can be declared multiple times corresponding to multiple '
        'control samples; see examples below'
    )
    samp_args.add_argument(
        '-x', '--ctrl-max', metavar='X', type=int, default=1,
        help='k-mers with abund > X in any control sample are uninteresting; '
        'default is X=1'
    )

    mem_desc = """\
    Specify how much memory to allocate for the sketch data structures
    used to store k-mer counts. The first control sample will be
    allocated the full amount of specifed `--memory`, and all subsequent
    samples will be allocated a fraction thereof.
    """
    mem_desc = textwrap.dedent(mem_desc)
    memory_args = subparser.add_argument_group('Memory allocation', mem_desc)
    memory_args.add_argument(
        '-M', '--memory', default='1e6', metavar='MEM',
        type=khmer_args.memory_setting, help='total memory to allocate for '
        'the initial control sample; default is 1M'
    )
    memory_args.add_argument(
        '-f', '--mem-frac', type=float, default=0.1, metavar='F',
        help='fraction of the total memory to allocate to subsequent samples; '
        'default is 0.1'
    )
    memory_args.add_argument(
        '--max-fpr', type=float, default=0.2, metavar='FPR',
        help='terminate if the expected false positive rate for any sample is '
        'higher than the specified FPR; default is 0.2'
    )

    band_desc = """\
    If memory is a limiting factor, it is possible to get a linear
    decrease in memory consumption by running kevlar in "banded" mode.
    Splitting the hashed k-mer space into N bands and only considering k-mers
    from one band at a time reduces the memory consumption to approximately 1/N
    of the total memory required. This implements a scatter/gather approach in
    which `kevlar count` and/or `kevlar novel` is run N times, after which the
    results are combined using `kevlar filter`.
    """
    band_desc = textwrap.dedent(band_desc)
    band_args = subparser.add_argument_group('K-mer banding', band_desc)
    band_args.add_argument(
        '--num-bands', type=int, metavar='N', default=None,
        help='number of bands into which to divide the hashed k-mer space'
    )
    band_args.add_argument(
        '--band', type=int, metavar='I', default=None,
        help='a number between 1 and N (inclusive) indicating the band to be '
        'processed'
    )

    misc_args = subparser.add_argument_group('Miscellaneous settings')
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
