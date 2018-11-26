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
    """Define the `kevlar novel` command-line interface."""

    desc = """\
    Identify "interesting" (potentially novel) k-mers and output the
    corresponding reads. Here we define "interesting" k-mers as those which are
    high abundance in each case sample and effectively absent (below some
    specified abundance threshold) in each control sample."""
    desc = textwrap.dedent(desc)
    epilog = """\
    Example::

        kevlar novel --out novel-reads.augfastq --case proband-reads.fq.gz \\
            --control father-reads-r1.fq.gz father-reads-r2.fq.gz \\
            --control mother-reads.fq.gz

    Example::

        kevlar novel --out novel-reads.augfastq.gz \\
            --control-counts father.counttable mother.counttable \\
            --case-counts proband.counttable --case proband-reads.fastq \\
            --ctrl-max 0 --case-min 10 --ksize 27

    Example::

        kevlar novel --out output.augfastq \\
            --case proband1.fq --case proband2.fq \\
            --control control1a.fq control1b.fq \\
            --control control2a.fq control2b.fq \\
            --save-case-counts p1.ct p2.ct --save-ctrl-counts c1.ct c2.ct"""
    epilog = textwrap.dedent(epilog)
    subparser = subparsers.add_parser(
        'novel', description=desc, epilog=epilog, add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    samp_desc = """\
    Specify input files, as well as thresholds for selecting "interesting"
    k-mers. A single pass is made over input files for control samples (to
    compute k-mer abundances), while two passes are made over input files for
    case samples (to compute k-mer abundances, and then to identify
    "interesting" k-mers). The k-mer abundance computing steps can be skipped
    if pre-computed k-mer abunandances are provided using the "--case-counts"
    and/or "--control-counts" settings. If "--control-counts" is declared, then
    all "--control" flags are ignored. If "--case-counts" is declared,
    FASTA/FASTQ files must still be provided with "--case" for selecting
    "interesting" k-mers and reads."""
    samp_desc = textwrap.dedent(samp_desc)
    samp_args = subparser.add_argument_group('Case/control config', samp_desc)
    samp_args.add_argument(
        '--case', metavar='F', nargs='+', required=True, action='append',
        help='one or more FASTA/FASTQ files containing reads from a case '
        'sample; can be declared multiple times corresponding to multiple '
        'case samples, see examples below'
    )
    samp_args.add_argument(
        '--case-counts', metavar='F', nargs='+',
        help='counttable file(s) corresponding to each case sample; if not '
        'provided, k-mer abundances will be computed from FASTA/FASTQ input; '
        'only one counttable per sample, see examples below'
    )
    samp_args.add_argument(
        '--control', metavar='F', nargs='+', action='append',
        help='one or more FASTA/FASTQ files containing reads from a control '
        'sample; can be declared multiple times corresponding to multiple '
        'control samples, see examples below'
    )
    samp_args.add_argument(
        '--control-counts', metavar='F', nargs='+',
        help='counttable file(s) corresponding to each control sample; if not '
        'provided, k-mer abundances will be computed from FASTA/FASTQ input; '
        'only one counttable per sample, see examples below'
    )
    samp_args.add_argument(
        '-x', '--ctrl-max', metavar='X', type=int, default=1,
        help='k-mers with abund > X in any control sample are uninteresting; '
        'default is X=1'
    )
    samp_args.add_argument(
        '-y', '--case-min', metavar='Y', type=int, default=6,
        help='k-mers with abund < Y in any case sample are uninteresting; '
        'default is Y=6'
    )
    samp_args.add_argument(
        '-M', '--memory', default='1e6', type=khmer_args.memory_setting,
        metavar='MEM',
        help='total memory allocated to k-mer abundance for each sample; '
        'default is 1M; ignored when pre-computed k-mer abundances are '
        'supplied via counttable'
    )
    samp_args.add_argument(
        '--max-fpr', type=float, default=0.2, metavar='FPR',
        help='terminate if the expected false positive rate for any sample is '
        'higher than the specified FPR; default is 0.2'
    )

    band_desc = """\
    If memory is a limiting factor, it is possible to get a linear decrease in
    memory consumption by running `kevlar novel` in "banded" mode. Splitting
    the hashed k-mer space into N bands and only considering k-mers from one
    band at a time reduces the memory consumption to approximately 1/N of the
    total memory required. This implements a scatter/gather approach in which
    `kevlar novel` is run N times, after the results are combined using
    `kevlar filter`."""
    band_desc = textwrap.dedent(band_desc)
    band_args = subparser.add_argument_group('K-mer banding', band_desc)
    band_args.add_argument('--num-bands', type=int, metavar='N', default=None,
                           help='number of bands into which to divide the '
                           'hashed k-mer space')
    band_args.add_argument('--band', type=int, metavar='I', default=None,
                           help='a number between 1 and N (inclusive) '
                           'indicating the band to be processed')

    out_args = subparser.add_argument_group('Output settings')
    out_args.add_argument('-o', '--out', metavar='FILE',
                          help='file to which interesting reads will be '
                          'written; default is terminal (stdout)')
    out_args.add_argument('--save-case-counts', metavar='CT', nargs='+',
                          help='save the computed k-mer counts for each case '
                          'sample to the specified count table file(s)')
    out_args.add_argument('--save-ctrl-counts', metavar='CT', nargs='+',
                          help='save the computed k-mer counts for each '
                          'control sample to the specified count table '
                          'file(s)')

    misc_args = subparser.add_argument_group('Miscellaneous settings')
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    misc_args.add_argument('--abund-screen', type=int, default=None,
                           metavar='INT', help='discard reads with any k-mers '
                           'whose abundance is < INT')
    misc_args.add_argument('-t', '--threads', type=int, default=1, metavar='T',
                           help='number of threads to use for file processing;'
                           ' default is 1')
    misc_args.add_argument('--skip-until', type=str, metavar='ID',
                           help='when re-running `kevlar novel`, skip all '
                           'reads in the case input until read with name `ID` '
                           'is observed')
