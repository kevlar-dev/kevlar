#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import kevlar
import screed
from kevlar import KmerOfInterest


def test_calc_offset_same_orientation():
    """
    GCTGCACCGATGTACGCAAA
                  |||||
                 ACGCAAAGCTATTTAAAACC
    """
    record1 = screed.Record(
        name='read1',
        sequence='GCTGCACCGATGTACGCAAA',
        ikmers=[KmerOfInterest('CGCAA', 14, [15, 0, 0])],
    )
    record2 = screed.Record(
        name='read2',
        sequence='ACGCAAAGCTATTTAAAACC',
        ikmers=[KmerOfInterest('CGCAA', 1, [15, 0, 0])],
    )
    tail, head, offset, sameorient = kevlar.assemble.calc_offset(
        record1, record2, 'CGCAA'
    )
    assert tail == record1
    assert head == record2
    assert offset == 13
    assert sameorient is True


def test_calc_offset_opposite_orientation():
    """
    GCTGCACCGATGTACGCAAA
                  |||||
                 ACGCAAAGCTATTTAAAACC <-- reverse complement
    """
    record1 = screed.Record(
        name='read1',
        sequence='GCTGCACCGATGTACGCAAA',
        ikmers=[KmerOfInterest('CGCAA', 14, [15, 0, 0])],
    )
    record2 = screed.Record(
        name='read2',
        sequence='GGTTTTAAATAGCTTTGCGT',
        ikmers=[KmerOfInterest('TTGCG', 14, [15, 0, 0])],
    )
    tail, head, offset, sameorient = kevlar.assemble.calc_offset(
        record1, record2, 'CGCAA'
    )
    assert tail == record1
    assert head == record2
    assert offset == 13
    assert sameorient is False
