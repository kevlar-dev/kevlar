#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import filecmp
import json
import sys
from tempfile import NamedTemporaryFile

import pandas
import pytest
import kevlar
from kevlar.dist import count_first_pass, count_second_pass
from kevlar.dist import calc_mu_sigma, compute_dist
from kevlar.dist import KevlarZeroAbundanceDistError
from kevlar.tests import data_file
from khmer import Counttable, Nodetable


def test_count_first_pass():
    mask = Nodetable.load(data_file('minitrio/mask.nt'))
    counts = Counttable(31, 1e4, 4)
    seqfile = data_file('minitrio/trio-proband.fq.gz')
    count_first_pass([seqfile], counts, mask)
    with NamedTemporaryFile(suffix='.ct') as countfile:
        counts.save(countfile.name)
        testcountfile = data_file('minitrio/trio-proband-mask-counts.ct')
        assert filecmp.cmp(testcountfile, countfile.name) is True


def test_count_second_pass():
    mask = Nodetable.load(data_file('minitrio/mask.nt'))
    counts = Counttable.load(data_file('minitrio/trio-proband-mask-counts.ct'))
    seqfile = data_file('minitrio/trio-proband.fq.gz')
    abund = count_second_pass([seqfile], counts)
    assert abund == {10: 6, 11: 10, 12: 12, 13: 18, 14: 16, 15: 11, 16: 9,
                     17: 9, 18: 11, 19: 8, 20: 9, 21: 7, 22: 3}


def test_calc_mu_sigma():
    abund = {10: 6, 11: 10, 12: 12, 13: 18, 14: 16, 15: 11, 16: 9,
             17: 9, 18: 11, 19: 8, 20: 9, 21: 7, 22: 3}
    mu, sigma = calc_mu_sigma(abund)
    assert pytest.approx(15.32558, mu)
    assert pytest.approx(3.280581, sigma)


def test_musigma_empty_dist():
    with pytest.raises(KevlarZeroAbundanceDistError):
        empty = dict()
        calc_mu_sigma(empty)


def test_compute_dist():
    abund = {10: 6, 11: 10, 12: 12, 13: 18, 14: 16, 15: 11, 16: 9,
             17: 9, 18: 11, 19: 8, 20: 9, 21: 7, 22: 3}
    data = compute_dist(abund)
    print(data)
    assert list(data['Count'][:5]) == [6.0, 10.0, 12.0, 18.0, 16.0]
    assert list(data['CumulativeCount'][:5]) == [6.0, 16.0, 28.0, 46.0, 62.0]


def test_dist():
    mask = Nodetable.load(data_file('minitrio/mask.nt'))
    filenames = [data_file('minitrio/trio-proband.fq.gz')]
    mu, sigma, data = kevlar.dist.dist(filenames, mask, memory=4e4)
    assert pytest.approx(15.32558, mu)
    assert pytest.approx(3.280581, sigma)
    assert list(data['Count'][-5:]) == [11.0, 8.0, 9.0, 7.0, 3.0]


def test_dist_empty():
    mask = Nodetable(31, 1e4, 4)
    mask.consume('GATTACA' * 10)
    mask.consume('A' * 50)
    filenames = [data_file('minitrio/trio-proband.fq.gz')]
    with pytest.raises(KevlarZeroAbundanceDistError):
        kevlar.dist.dist(filenames, mask, memory=4e4)


def test_main(capsys):
    arglist = [
        'dist', data_file('minitrio/mask.nt'),
        data_file('minitrio/trio-proband.fq.gz')
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.dist.main(args)

    out, err = capsys.readouterr()
    print(out)
    js = json.loads(out)
    assert pytest.approx(15.32558, js['mu'])
    assert pytest.approx(3.280581, js['sigma'])


def test_main_plot():
    with NamedTemporaryFile(suffix='.png') as plotfile:
        arglist = [
            'dist', '--plot', plotfile.name, '--plot-xlim', '0', '50',
            data_file('minitrio/mask.nt'),
            data_file('minitrio/trio-proband.fq.gz')
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.dist.main(args)


def test_tsv():
    with NamedTemporaryFile(suffix='.tsv') as distfile:
        arglist = [
            'dist', '--tsv', distfile.name, data_file('minitrio/mask.nt'),
            data_file('minitrio/trio-proband.fq.gz')
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.dist.main(args)
        data = pandas.read_table(distfile.name, sep='\t')
        print(data)
        cuml_test = [
            15.0, 18.0, 24.0, 44.0, 78.0, 153.0, 222.0, 325.0, 423.0, 515.0,
            585.0, 666.0, 756.0, 814.0, 861.0, 888.0, 902.0, 903.0
        ]
        assert list(data['CumulativeCount']) == cuml_test
