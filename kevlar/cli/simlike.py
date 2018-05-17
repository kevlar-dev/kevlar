#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys


def subparser(subparsers):
    """Define the `kevlar simlike` command-line interface."""

    desc = "Sort variants by likelihood score"
    subparser = subparsers.add_parser('simlike', description=desc)
    subparser.add_argument('--case', metavar='CT',
                           help='k-mer counttable for case/proband')
    subparser.add_argument('--controls', nargs='+', metavar='CT',
                           help='k-mer counttables for controls/parents/'
                           'siblings, 1 counttable per sample')
    subparser.add_argument('--mu', metavar='μ', type=float, default=30.0,
                           help='Mean k-mer abundance; default is 30.0')
    subparser.add_argument('--sigma', metavar='σ', type=float, default=8.0,
                           help='Standard deviation of k-mer abundance; '
                           'default is 8.0')
    subparser.add_argument('--epsilon', metavar='ε', type=float, default=0.001,
                           help='Error rate; default is 0.001')
    subparser.add_argument('-o', '--out', metavar='OUT', default=sys.stdout,
                           help='output file; default is terminal (standard '
                           'output)')
    subparser.add_argument('--case-min', metavar='Y', type=int, default=5,
                           help='Minimum abundance threshold for proband; '
                           'default is 5')
    subparser.add_argument('--sample-labels', metavar='LBL', type=str,
                           nargs='+', help='List of sample labels (with '
                           'proband first)')
    subparser.add_argument('--refr', metavar='REFR', help='Small counttable of'
                           ' k-mer counts in the reference genome')
    subparser.add_argument('vcf')
