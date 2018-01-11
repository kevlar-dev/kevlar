#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    """Define the `kevlar gentrio` command-line interface."""

    desc = 'Apply randomly generated mutations to the genome provided.'

    subparser = subparsers.add_parser('gentrio', description=desc)
    subparser.add_argument('-i', '--inherited', type=int, metavar='I',
                           default=20, help='number of shared/inherited '
                           'mutations to simulate')
    subparser.add_argument('-d', '--de-novo', type=int, metavar='D',
                           default=10, help='number of unique/de novo '
                           'mutations to simulate')
    subparser.add_argument('--vcf', metavar='FILE', help='write mutations to '
                           'a VCF file')
    subparser.add_argument('--prefix', metavar='PFX', default='trio',
                           help='prefix for output fasta files; default is '
                           '"trio"')
    subparser.add_argument('--weights', metavar='WT',
                           default='snv=0.8,ins=0.1,del=0.1',
                           help='comma-separated list of key/value pairs '
                           'indicating the relative frequency of different '
                           'variant types; default is "snv=0.8,ins=0.1,'
                           'del=0.1"; ')
    subparser.add_argument('-s', '--seed', metavar='S', default=None, type=int,
                           help='seed for random number generator')
    subparser.add_argument('genome', help='genome to mutate')
