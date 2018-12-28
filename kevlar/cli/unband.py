#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 Battelle National Biodefense Institute
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
import textwrap
import khmer
from khmer import khmer_args


def subparser(subparsers):
    """Define the `kevlar unband` command-line interface."""

    desc = """\
    When kevlar is run in k-mer banding mode, the same read will typically
    appear in multiple output files, annotated with a different set of
    potentially novel k-mers in each case. This command will consolidate these
    duplicated records across files into a single non-redundant set of reads
    with the complete set of novel k-mer annotations.
    """
    desc = textwrap.dedent(desc)

    subparser = subparsers.add_parser(
        'unband', description=desc,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparser.add_argument(
        '-n', '--n-batches', metavar='N', type=int, default=16,
        help='number of batches into which records will be split; default is '
        '16; N temporary files are created and each record from the input is '
        'written to a temporary file based on its read name; then each batch '
        'is loaded into memory and duplicate records are resolved'
    )
    subparser.add_argument(
        '-o', '--out', metavar='FILE',
        help='output file; default is terminal (stdout)'
    )
    subparser.add_argument(
        'infile', nargs='+', help='input files in augmented Fasta/Fastq format'
    )
