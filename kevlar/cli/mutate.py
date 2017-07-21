#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    """Define the `kevlar mutate` command-line interface."""

    desc = 'Apply the specified mutations to the genome provided.'

    subparser = subparsers.add_parser('mutate', description=desc)
    subparser.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('mutations', help='mutations file')
    subparser.add_argument('genome', help='genome to mutate')
