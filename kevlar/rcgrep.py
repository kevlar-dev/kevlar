#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import argparse
import gzip
import sys

import kevlar


def subparser(subparsers):
    subparser = subparsers.add_parser('rcgrep')
    subparser.add_argument('-q', '--query', action='append', metavar='SEQ',
                           required=True)
    subparser.add_argument('file', nargs='+')


def main(args):
    queryseqs = set()
    for qseq in args.query:
        queryseqs.add(qseq)
        queryseqs.add(kevlar.revcom(qseq))

    for infile in args.file:
        openfunc = gzip.open if infile.endswith('.gz') else open
        instream = openfunc(infile, 'rt')
        for line in instream:
            for seq in queryseqs:
                if seq in line:
                    print(line, end='')
                    break
