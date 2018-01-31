#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import sys
import khmer
import kevlar
from math import log
from scipy.stats import norm
from kevlar.tests import data_file
from kevlar.simlike import filter_refr, set_error_rates, abund_log_prob


def test_filter_refr():
    """
    refr     AATTATTGTTGTATGTATGTGCAAACTCTTAATTATACTAGCAATATATTT
    alt      AATTATTGTTGTATGTATGTGCAA---CTTAATTATACTAGCAATATATTT
    """
    refrcounts = khmer.SmallCounttable(25, 1e4, 4)
    refrseq = 'AATTATTGTTGTATGTATGTGCAAACTCTTAATTATACTAGCAATATATTT'
    altseq = 'AATTATTGTTGTATGTATGTGCAACTTAATTATACTAGCAATATATTT'
    altkmers = refrcounts.get_kmers(altseq)
    refrkmers = refrcounts.get_kmers(refrseq)
    newalt, newrefr = kevlar.simlike.filter_refr(altkmers, refrkmers,
                                                 refrcounts)
    assert altkmers == newalt
    assert refrkmers == newrefr

    r = 'AATCTAAGAAATTATTGTTGTATGTATGTGCAAACTCTTAATTATACTAGCAATATATTTGAGAGATGG'
    refrcounts.consume(r)
    newalt, newrefr = kevlar.simlike.filter_refr(altkmers, refrkmers,
                                                 refrcounts)
    assert altkmers == newalt
    assert refrkmers == newrefr

    refrcounts.consume('TTGTTGTATGTATGTGCAAACTCTTAATT')
    newalt, newrefr = kevlar.simlike.filter_refr(altkmers, refrkmers,
                                                 refrcounts)
    assert altkmers == newalt
    assert len(refrkmers) == len(newrefr) + 5


def test_get_abund():
    casecounts = khmer.Counttable(25, 1e4, 4)
    ctrl1counts = khmer.Counttable(25, 1e4, 4)
    ctrl2counts = khmer.Counttable(25, 1e4, 4)

    dna = 'AATTATTGTTGTATGTATGTGCAAACTCTTAATTATACTAGCAATATATTT'

    for _ in range(5):
        casecounts.consume(dna)
    for _ in range(2):
        ctrl1counts.consume(dna)
    for _ in range(15):
        ctrl2counts.consume(dna)

    kmers = casecounts.get_kmers(dna)
    controls = (ctrl1counts, ctrl2counts)
    abunds = kevlar.simlike.get_abundances(kmers, casecounts, controls)
    print(*abunds, sep='\n')
    testabunds = [
        [5] * len(kmers),
        [2] * len(kmers),
        [15] * len(kmers),
    ]
    assert abunds == testabunds


def test_set_error_rates():
    assert set_error_rates(0.01, 3) == [0.01, 0.01, 0.01]
    assert set_error_rates([0.01, 0.02, 0.03], 3) == [0.01, 0.02, 0.03]
    with pytest.raises(ValueError) as ve:
        set_error_rates(None, 4)
    assert "doesn't look/smell/quack like a float or list of floats" in str(ve)


def test_abund_log_prob():
    for abund in (1, 3, 47, 255):
        assert abund_log_prob(0, abund, error=0.01) == abund * log(0.01)
    assert abund_log_prob(1, 30, mean=30.0, sd=6.0) == norm.logpdf(30, 15, 3)
    assert abund_log_prob(2, 30, mean=40.0, sd=8.0) == norm.logpdf(30, 40, 8)
    with pytest.raises(kevlar.simlike.KevlarUnknownGenotypeError) as ge:
        abund_log_prob('SUPDOG', 15)
