#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import sys
import kevlar
from kevlar.tests import data_file


@pytest.mark.parametrize('greedy', [True, False])
def test_pico_4(greedy, capsys):
    reads = data_file('pico-4.augfastq.gz')
    refr = data_file('human-random-pico.fa.gz')
    arglist = ['alac', '--ksize', '25', reads, refr]
    args = kevlar.cli.parser().parse_args(arglist)
    args.greedy = greedy
    kevlar.alac.main(args)
    out, err = capsys.readouterr()

    vcf = '\t'.join([
        'seq1', '1175768', '.', 'T', 'C', '.', 'PASS',
        'RW=CCCTGCCATTATAGATGCTAGATTTACATCTTCATTTATTTTTACTTTT;'
        'VW=CCCTGCCATTATAGATGCTAGATTCACATCTTCATTTATTTTTACTTTT'
    ])
    assert vcf.strip() == out.strip()


@pytest.mark.parametrize('cc,pos,ref,alt', [
    (2, 834645, 'A', 'AGTGGGATTACGTAGGAAATCCGCGGGGCTGTGACATATATTTGTTGACAAGCATA'
                     'TATTGTTCCTAGAGGTCGTTGGGTTCGTTACACCCAAGGGGGCGTATAACATGTTA'
                     'CTCAGTTGCGTCGGACCGATTAATAACTCGAATGTAAGGCAGGATATTT'),
    (3, 4072, 'G', 'GCCGAGACGCAGCGTGATACTTAAGATTAAGTTAAGCAACAGCTTAGCGTACGCAATT'
                   'GCGTCTAATTGAGGGGCCGTAGATATAAGCTCCGTGTTCTCAGTTGGTGGGTAACAGA'
                   'ACCCGCAAGCACACCGCTTTCAGTGTGTCACATGCACA'),
    (5, 1175767, 'T', 'C'),
    (6, 185751, 'TCAAACTCTGGCATTATACATAGGGTTCCCG', 'T'),
    (7, 2265794, 'GCAGGGTACATAAGAGTCCATTGTGCCTGTATTATTTTGAGCAATGGCTAAAGTACCTTC'
                 'ACCCTTGCTC', 'G'),
    (8, 636698, 'C', 'A'),
    (9, 226610, 'TTCAACTCTACAGGGTCTGATGCTTACAGGAGTTCCCTTTTCCTACATTTGGTTCAAGATG'
                'GCAACAAATACATTTTAGATTCACATAGCTCATCCTTCTAGGTTAACAGTAAACTTAAGAA'
                'CTAAGACCAGAACCAGGAGGGTCAGGAAATTCTCCTGTGTGGTTGCTGGGACCACTGCAAA'
                'GCAGTGGC', 'T'),
    (10, 1527138, 'C', 'CTCCTGGTCTGCCACGGTTGACTTGCCTACATAT'),
])
def test_pico_calls(cc, pos, ref, alt):
    reads = data_file('pico-var/cc{:d}.afq.gz'.format(cc))
    readstream = kevlar.parse_augmented_fastx(kevlar.open(reads, 'r'))
    refrfile = data_file('human-random-pico.fa.gz')
    caller = kevlar.alac.alac(readstream, refrfile, ksize=25, delta=50)
    calls = [v for v in caller]

    assert len(calls) == 1
    assert calls[0]._pos == pos
    assert calls[0]._refr == ref
    assert calls[0]._alt == alt
