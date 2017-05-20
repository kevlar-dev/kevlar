#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import screed
import kevlar
from kevlar import KmerOfInterest
from kevlar.overlap import (calc_offset, OverlappingReadPair,
                            INCOMPATIBLE_PAIR)


@pytest.fixture
def record1():
    return screed.Record(
        name='read1',
        sequence='GCTGCACCGATGTACGCAAA',
        ikmers=[KmerOfInterest('CGCAA', 14, [15, 0, 0])],
    )


@pytest.fixture
def record2():
    return screed.Record(
        name='read2',
        sequence='ACGCAAAGCTATTTAAAACC',
        ikmers=[
            KmerOfInterest('CGCAA', 1, [15, 0, 0]),
            KmerOfInterest('AAAAC', 14, [19, 1, 0]),
        ],
    )


@pytest.fixture
def record3():
    # reverse complement of record2
    return screed.Record(
        name='read3',
        sequence='GGTTTTAAATAGCTTTGCGT',
        ikmers=[
            KmerOfInterest('GTTTT', 1, [19, 1, 0]),
            KmerOfInterest('TTGCG', 14, [15, 0, 0]),
        ],
    )


@pytest.fixture
def record4():
    # similar to record2 but with a single nucleotide mismatch
    return screed.Record(
        name='read4',
        sequence='ACGCAATGCTATTTAAAACC',
        ikmers=[
            KmerOfInterest('CGCAA', 1, [15, 0, 0]),
            KmerOfInterest('AAAAC', 14, [19, 1, 0]),
        ],
    )


@pytest.fixture
def record5():
    return screed.Record(
        name='read5',
        sequence='CTCTTCCGGCAGTCACTGTCAAGAGAGGGTGAACT',
        ikmers=[
            KmerOfInterest('CTGTCAA', 15, [12, 0, 0]),
            KmerOfInterest('TGTCAAG', 16, [13, 0, 0]),
        ],
    )


@pytest.fixture
def record6():
    return screed.Record(
        name='read6',
        sequence='TCACTGTCAAGAGAGGCCTACGGATTCGGTTACTG',
        ikmers=[
            KmerOfInterest('CTGTCAA', 3, [12, 0, 0]),
            KmerOfInterest('TGTCAAG', 4, [13, 0, 0]),
        ],
    )


def test_calc_offset_same_orientation(record1, record2):
    """
    Compute offset of reads sharing an interesting k-mer, same orientation.

    GCTGCACCGATGTACGCAAA
                  |||||
                 ACGCAAAGCTATTTAAAACC
    """
    pair = kevlar.overlap.calc_offset(record1, record2, 'CGCAA')
    assert pair.tail == record1
    assert pair.head == record2
    assert pair.offset == 13
    assert pair.overlap == 7
    assert pair.sameorient is True


def test_calc_offset_opposite_orientation(record1, record3):
    """
    Compute offset of reads sharing an interesting k-mer, opposite orientation.

    GCTGCACCGATGTACGCAAA
                  |||||
                 ACGCAAAGCTATTTAAAACC <-- reverse complement
    """
    pair = kevlar.overlap.calc_offset(record1, record3, 'CGCAA')
    assert pair.tail == record1
    assert pair.head == record3
    assert pair.offset == 13
    assert pair.overlap == 7
    assert pair.sameorient is False


def test_calc_offset_mismatch(record1, record4):
    """
    Compute offset of reads sharing an interesting k-mer, but with mismatch.

    GCTGCACCGATGTACGCAAA
                  |||||X
                 ACGCAATGCTATTTAAAACC

    The interesting k-mer is simply a seed. The assembler requires that the
    entire overlap between the two reads matches exactly, so the mismatch here
    should return a null offset indicating that the pair of reads is
    incompatible despite sharing an interesting k-mer.
    """
    pair = kevlar.overlap.calc_offset(record1, record4, 'CGCAA')
    assert pair == kevlar.overlap.INCOMPATIBLE_PAIR


def test_calc_offset_weirdness(record5, record6):
    """
    Compute offset of reads sharing an interesting k-mer, but with mismatches.

    CTCTTCCGGCAGTCACTGTCAAGAGAGGGTGAACT
                   |||||||
                    |||||||     XXXXXXX
                TCACTGTCAAGAGAGGCCTACGGATTCGGTTACTG

    Not just a single mismatch, but extensive differences making the reads
    incompatible for assembly.
    """
    for ikmer in ['CTGTCAA', 'TGTCAAG']:
        pair = kevlar.overlap.calc_offset(record5, record6, ikmer)
        assert pair == kevlar.overlap.INCOMPATIBLE_PAIR
