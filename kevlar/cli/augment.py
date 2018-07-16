#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    """Define the `kevlar augment` command-line interface."""

    desc = """\
    Internally, kevlar annotates sequences with "interesting k-mers" and uses
    "augmented" Fastq and Fasta formats. Processing sequences with third-part
    tools usually requires discarding these annotations. This command is used
    to augment/reaugment a set of sequences using annotations from an already
    augmented sequence file.
    """

    subparser = subparsers.add_parser('augment', description=desc)
    subparser.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('-c', '--collapse-mates', action='store_true',
                           help='annotate output sequence with all mate '
                           'seqeuences')
    subparser.add_argument('augseqs', help='augmented sequence file')
    subparser.add_argument('seqs', help='sequences to annotate')
