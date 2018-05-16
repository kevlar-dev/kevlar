#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import khmer
import kevlar
from kevlar.tests import data_file
from kevlar.simlike import get_abundances, abund_log_prob
import pytest


def test_get_abundances():
    kid = khmer.Counttable(31, 1e6, 4)
    mom = khmer.Counttable(31, 1e6, 4)
    dad = khmer.Counttable(31, 1e6, 4)
    ref = khmer.Counttable(31, 1e6, 4)
    kid.consume_seqfile(data_file('minitrio/trio-proband.fq.gz'))
    mom.consume_seqfile(data_file('minitrio/trio-mother.fq.gz'))
    dad.consume_seqfile(data_file('minitrio/trio-father.fq.gz'))
    ref.consume_seqfile(data_file('minitrio/refr.fa'))

    window = 'ACTACCCTAACTTTGGAAGACTATCAAAAACCCATTTCTGGGGGTGGAGGGGAGGGAGACA'
    abund, ndropped = get_abundances(window, kid, (mom, dad), ref)
    assert ndropped == 0
    assert abund == [
        [7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 7, 7, 9, 9, 9, 9, 8, 8, 9, 7, 6, 6, 6,
         6, 6, 6, 6, 6, 6, 6, 7],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
         1, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1, 0, 0, 0, 0, 0, 0, 0],
    ]


def test_abund_log_prob():
    assert abund_log_prob(0, 3, refrabund=1) == pytest.approx(-10.51967)
    assert abund_log_prob(0, 4, refrabund=1) == pytest.approx(-14.02623)
    assert abund_log_prob(0, 4, refrabund=6) == pytest.approx(-6.85919)
    assert abund_log_prob(0, 4, refrabund=15) == pytest.approx(-3.19403)

    assert abund_log_prob(1, 1) == pytest.approx(-8.43023)
    assert abund_log_prob(1, 10) == pytest.approx(-3.08648)
    assert abund_log_prob(1, 15) == pytest.approx(-2.305233)
    assert abund_log_prob(1, 20) == pytest.approx(-3.08648)
    assert abund_log_prob(1, 10, mean=50.0, sd=9.9) == pytest.approx(-7.10969)
    assert abund_log_prob(1, 20, mean=50.0, sd=9.9) == pytest.approx(-3.02848)

    assert abund_log_prob(2, 1) == pytest.approx(-9.5687)
    assert abund_log_prob(2, 10) == pytest.approx(-6.12338)
    assert abund_log_prob(2, 30) == pytest.approx(-2.99838)
    assert abund_log_prob(2, 53) == pytest.approx(-7.13119)
    assert abund_log_prob(2, 29, mean=47.0, sd=9.3) == pytest.approx(-5.0220)
    assert abund_log_prob(2, 37, mean=47.0, sd=9.3) == pytest.approx(-3.727054)
    assert abund_log_prob(2, 43, mean=47.0, sd=9.3) == pytest.approx(-3.241449)
