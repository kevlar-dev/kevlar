#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
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
    Compute k-mer abundances for the provided samples.

    The `effcount` subcommand uses a memory-efficient strategy in which the
    first sample loaded acts as a mask for all subsequent samples. Any k-mer
    present in the mask (aboved some threshold) is not stored in any subsequent
    samples. As a result, many uninteresting k-mers can be ignored and the
    potentially interesting k-mers stored accurately in much less memory.

    Also supports k-mer banding. See
    http://kevlar.readthedocs.io/en/latest/banding.html for more details
    """
    desc = textwrap.dedent(desc)

    epilog = """\
    Example::

        kevlar effcount \\
            --sample control-reads-R1.fastq control-reads-R2.fastq \\
            --sample case-reads-R1.fastq case-reads-R2.fastq \\
            --memory 500M \\
            case1.ct case2.ct

    Example::

        kevlar effcount \\
            --sample mother-reads.fq.gz \\
            --sample father-reads.fq.gz \\
            --sample proband-R1.fq.gz proband-R2.fq.gz \\
            --ksize 25 --memory 12G --memfrac 0.1 --max-fpr 0.01 --threads 8 \\
            proband.counttable mother.counttable father.counttable"""
    epilog = textwrap.dedent(epilog)

    subparser = subparsers.add_parser(
        'effcount', description=desc, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparser.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    subparser.add_argument('-M', '--memory', type=khmer_args.memory_setting,
                           default=1e6, metavar='MEM',
                           help='memory to allocate for the initial count '
                           'table')
    subparser.add_argument('--memfrac', type=float, default=0.25, metavar='F',
                           help='proportion of `MEM` to allocate to all '
                           'subsequent samples; default is 0.25')
    subparser.add_argument('--max-fpr', type=float, default=0.2,
                           metavar='FPR', help='terminate if the estimated '
                           'false positive rate for any sample is higher than '
                           '"FPR"; default is 0.2')
    subparser.add_argument('-x', '--max-abund', type=int, default=5,
                           metavar='A', help='k-mers with abundance >= A in '
                           'the first sample are ignored in all subsequent '
                           'samples; default is 1')
    subparser.add_argument('--num-bands', type=int, metavar='N', default=None,
                           help='number of bands into which to divide the '
                           'hashed k-mer space')
    subparser.add_argument('--band', type=int, metavar='I', default=None,
                           help='a number between 1 and N (inclusive) '
                           'indicating the band to be processed')
    subparser.add_argument('-t', '--threads', type=int, default=1, metavar='T',
                           help='number of threads to use for file processing;'
                           ' default is 1')
    subparser.add_argument('--sample', type=str, nargs='+', metavar='FQ',
                           action='append', help='list of Fastq/Fasta input '
                           'files, declared separately for each sample')
    subparser.add_argument('outfiles', type=str, nargs='+', help='name of the '
                           'file to which the output (a k-mer count table) '
                           'will be written; the suffix ".counttable" will be '
                           'applied if the provided file name does not end in '
                           '".ct" or ".counttable"')
