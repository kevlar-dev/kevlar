#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import deque
import argparse
import gzip
import sys

import kevlar


def linesets(instream, before=0, after=0):
    """
    Yield a line along with its before-and-after context.
    """
    if before + after == 0:
        for line in instream:
            yield deque([line])
    else:
        lines = deque()
        for line in instream:
            lines.append(line)
            if len(lines) == before + after + 1:
                yield lines
                lines.popleft()


def positive_int(value):
    ival = int(value)
    if ival < 1:
        message = '"{}" is not a positive integer'.format(value)
        raise argparse.ArgumentTypeError(message)
    return ival


def subparser(subparsers):
    subparser = subparsers.add_parser('rcgrep')
    subparser.add_argument('-A', '--after', metavar='N', type=positive_int,
                           default=0)
    subparser.add_argument('-B', '--before', metavar='N', type=positive_int,
                           default=0)
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
        for lineset in linesets(instream, args.before, args.after):
            queryline = lineset[args.before]
            for seq in queryseqs:
                if seq in queryline:
                    for line in lineset:
                        print(line, end='')
                    break
