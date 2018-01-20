#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import khmer
from khmer import khmer_args


def subparser(subparsers):
    """Define the `kevlar genexp` command-line interface."""

    desc = 'Generate data for a case/control experiment.'

    subparser = subparsers.add_parser('genexp', description=desc)
    subparser.add_argument('genome', help='genome to mutate')
