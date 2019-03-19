#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    subparser = subparsers.add_parser(
        'simlike', description='Sort variants by likelihood score',
        add_help=False
    )

    count_args = subparser.add_argument_group(
        'K-mer count files',
        'Likelihood scores are based on the abundance of alternate allele '
        'k-mers in each sample and on the abundance of reference allele '
        'k-mers in the reference genome.'
    )
    count_args.add_argument(
        '--case', metavar='CT', required=True,
        help='k-mer counttable for case/proband'
    )
    count_args.add_argument(
        '--controls', nargs='+', metavar='CT', required=True,
        help='k-mer counttables for controls/parents/siblings, 1 counttable '
        'per sample'
    )
    count_args.add_argument(
        '--refr', metavar='REFR',  required=True,
        help='k-mer smallcounttable for reference genome'
    )

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

    filt_args = subparser.add_argument_group(
        'Heuristic filters',
        'The following heuristic filters can improve accuracy when calling '
        'de novo variants but may require tuning for your particular data set.'
    )
    filt_args.add_argument(
        '--ctrl-abund-high', metavar='H', type=int, default=4,
        help='a variant call will be filtered out if either of the control '
        'samples has >H high abundance k-mers spanning the variant (where '
        'high abundance means > --ctrl-max); by default, H=4; set H<=0 to '
        'disable the filter'
    )
    filt_args.add_argument(
        '--case-abund-low', metavar='L', type=int, default=5,
        help='a variant call will be filtered out if the case sample has L or '
        'more consecutive low abundance k-mers spanning the variant (where '
        'low abundance means < --case-min); by default, L=5; set L<=0 to '
        'disable the filter'
    )
    filt_args.add_argument(
        '--min-like-score', metavar='S', type=float, default=0.0,
        help='filter out variant predictions with likelihood scores < S; by '
        'default, S = 0.0, but it\'s often possible to improve specificity '
        'without sacrificing sensitivity by raising S to, for example, 50.0'
    )
    filt_args.add_argument(
        '--drop-outliers', action='store_true',
        help='discard terminal variant-spanning k-mers with abunance much '
        'higher than average (representing k-mers that should be in the '
        'reference genome but are not); this will increase sensitivity, but '
        'will potentially introduce many false calls as well'
    )
    filt_args.add_argument(
        '--ambig-thresh', metavar='A', type=int, default=10,
        help='discard contigs that result in > A distinct, equally optimal '
        'variant calls; by default, A = 10; set A=0 to disable this filter'
    )

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
