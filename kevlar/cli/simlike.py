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
    subparser = subparsers.add_parser('simlike', description=desc,
                                      add_help=False)

    count_args = subparser.add_argument_group(
        'K-mer count files',
        'Likelihood scores are based on the abundance of alternate allele '
        'k-mers in each sample and on the abundance of reference allele '
        'k-mers in the reference genome.'
    )
    count_args.add_argument('--case', metavar='CT',
                            help='k-mer counttable for case/proband')
    count_args.add_argument('--controls', nargs='+', metavar='CT',
                            help='k-mer counttables for controls/parents/'
                            'siblings, 1 counttable per sample')
    count_args.add_argument('--refr', metavar='REFR', help='k-mer '
                            'smallcounttable for reference genome')

    thresh_args = subparser.add_argument_group(
        'K-mer count thresholds',
        'The thresholds originally used to detect novel k-mers are used at '
        'this stage to distinguish true variants from spurious predictions.'
    )
    thresh_args.add_argument('--ctrl-max', metavar='X', type=int, default=1,
                             help='maximum abundance threshold for controls; '
                             'default is 1')
    thresh_args.add_argument('--case-min', metavar='Y', type=int, default=6,
                             help='minimum abundance threshold for proband; '
                             'default is 6')

    cov_args = subparser.add_argument_group(
        'K-mer coverage',
        'Likelihood scores also depend on the estimated or observed '
        'distribution of k-mer abundances in each sample.'
    )
    cov_args.add_argument('--mu', metavar='μ', type=float, default=30.0,
                          help='mean k-mer abundance; default is 30.0')
    cov_args.add_argument('--sigma', metavar='σ', type=float, default=8.0,
                          help='standard deviation of k-mer abundance; '
                          'default is 8.0')
    cov_args.add_argument('--epsilon', metavar='ε', type=float, default=0.001,
                          help='error rate; default is 0.001')
    cov_args.add_argument('--static', dest='dynamic', action='store_false',
                          help='when computing likelihood scores for SNVs, '
                          'disable dynamic error rate (scaled by abundance of '
                          'reference allele k-mer in the reference genome)')

    misc_args = subparser.add_argument_group('Miscellaneous settings')
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('--sample-labels', metavar='LBL', type=str,
                           nargs='+', help='list of sample labels (with case/'
                           'proband first)')
    misc_args.add_argument('-f', '--fast-mode', action='store_true',
                           help='whenever possible, stop computations '
                           'prematurely for any putative variants that have '
                           'already been filtered out')
    misc_args.add_argument('-o', '--out', metavar='OUT', default='-',
                           help='output file; default is terminal (standard '
                           'output)')
    subparser.add_argument('vcf', nargs='+')
