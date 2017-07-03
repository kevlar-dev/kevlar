#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import screed
import sys
import kevlar
from kevlar import KmerOfInterest
from kevlar.overlap import (print_read_pair, calc_offset, OverlappingReadPair,
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
def record2a():
    return screed.Record(
        name='read2',
        sequence='ACGCAAAGCTATTTACGCAA',
        ikmers=[
            KmerOfInterest('CGCAA', 1, [15, 0, 0]),
            KmerOfInterest('CGCAA', 15, [15, 0, 0]),
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


@pytest.fixture
def picorecord1():
    return screed.Record(
        name='seq1_901350_901788_1:0:0_0:0:0_21ca1/2',
        sequence=('GTTTTTTTTTTGTTTCCCAAAGTAAGGCTGAGTGAACAATATTTTCTCATAGTTTTGAC'
                  'AAAAACAAAGGAATCCTTAGTTATTAAACTCGGGAGTTTGA'),
        ikmers=[
            KmerOfInterest('TTTTTTGTTTCCCAAAGTAAGGCTG', 5, [19, 0, 0]),
            KmerOfInterest('TTTTTGTTTCCCAAAGTAAGGCTGA', 6, [18, 1, 0]),
            KmerOfInterest('TTTTGTTTCCCAAAGTAAGGCTGAG', 7, [18, 1, 0]),
            KmerOfInterest('TTTGTTTCCCAAAGTAAGGCTGAGT', 8, [18, 0, 0]),
            KmerOfInterest('TTGTTTCCCAAAGTAAGGCTGAGTG', 9, [17, 0, 0]),
        ],
    )


@pytest.fixture
def picorecord2():
    return screed.Record(
        name='seq1_901428_901847_3:0:0_0:0:0_87d/1',
        sequence=('TTACATTTATTCGTTTGTGCAGGCTGAGACCTCACTTCCAACTGTAATCCAAAAGCTTA'
                  'GTTTTTTTTTTGTTTCCCAAAGTAAGGCTGAGTGAACAATA'),
        ikmers=[
            KmerOfInterest('TTTTTTGTTTCCCAAAGTAAGGCTG', 64, [19, 0, 0]),
            KmerOfInterest('TTTTTGTTTCCCAAAGTAAGGCTGA', 65, [18, 1, 0]),
            KmerOfInterest('TTTTGTTTCCCAAAGTAAGGCTGAG', 66, [18, 1, 0]),
            KmerOfInterest('TTTGTTTCCCAAAGTAAGGCTGAGT', 67, [18, 0, 0]),
            KmerOfInterest('TTGTTTCCCAAAGTAAGGCTGAGTG', 68, [17, 0, 0]),
        ],
    )


@pytest.fixture
def picorecord3():
    return screed.Record(
        name='seq1_901428_901847_3:0:0_0:0:0_87d/1',
        sequence=('TATTGTTCACTCAGCCTTACTTTGGGAAACAAAAAAAAAACTAAGCTTTTGGATTACAG'
                  'TTGGAAGTGAGGTCTCAGCCTGCACAAACGAATAAATGTAA'),
        ikmers=[
            KmerOfInterest('CAGCCTTACTTTGGGAAACAAAAAA', 11, [17, 0, 0]),
            KmerOfInterest('TCAGCCTTACTTTGGGAAACAAAAA', 10, [18, 0, 0]),
            KmerOfInterest('CTCAGCCTTACTTTGGGAAACAAAA', 9, [18, 1, 0]),
            KmerOfInterest('ACTCAGCCTTACTTTGGGAAACAAA', 8, [18, 1, 0]),
            KmerOfInterest('CACTCAGCCTTACTTTGGGAAACAA', 7, [19, 0, 0]),
        ],
    )


def test_print_read_pair_same_orient(record1, record2, capsys):
    pair = kevlar.overlap.calc_offset(record1, record2, 'CGCAA')
    print_read_pair(pair, 14, sys.stderr)
    out, err = capsys.readouterr()

    assert 'read1 --(overlap=7, offset=13, sameorient=True)--> read2' in err
    assert ('GCTGCACCGATGTACGCAAA\n'
            '              |||||\n'
            '             ACGCAAAGCTATTTAAAACC') in err


def test_print_read_pair_diff_orient(record1, record3, capsys):
    pair = kevlar.overlap.calc_offset(record1, record3, 'CGCAA')
    print_read_pair(pair, 14, sys.stderr)
    out, err = capsys.readouterr()

    assert 'read1 --(overlap=7, offset=13, sameorient=False)--> read3' in err
    assert ('GCTGCACCGATGTACGCAAA\n'
            '              |||||\n'
            '             ACGCAAAGCTATTTAAAACC') in err


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


def test_calc_offset_kmer_multi_freq(record1, record2a):
    """
    Read contains multiple copies of the same k-mer.

    Make sure kevlar doesn't try to compute offset of reads when one of them
    contains more than one copy of the interesting k-mer.

    GCTGCACCGATGTACGCAAA
                  |||||         |||||
                 ACGCAAAGCTATTTACGCAA
    """
    pair = kevlar.overlap.calc_offset(record1, record2a, 'CGCAA')
    assert pair == INCOMPATIBLE_PAIR


def test_calc_offset_opposite_orientation(record1, record3, capsys):
    """
    Compute offset of reads sharing an interesting k-mer, opposite orientation.

    GCTGCACCGATGTACGCAAA
                  |||||
                 ACGCAAAGCTATTTAAAACC <-- reverse complement
    """
    pair = kevlar.overlap.calc_offset(record1, record3, 'CGCAA', sys.stderr)
    out, err = capsys.readouterr()
    assert pair.tail == record1
    assert pair.head == record3
    assert pair.offset == 13
    assert pair.overlap == 7
    assert pair.sameorient is False
    assert '--(overlap=7, offset=13, sameorient=False)-->' in err


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
    assert pair == INCOMPATIBLE_PAIR


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
        assert pair == INCOMPATIBLE_PAIR


def test_pico_offset(picorecord1, picorecord2, picorecord3):
    pair = kevlar.overlap.calc_offset(picorecord1, picorecord2,
                                      'TTTTTTGTTTCCCAAAGTAAGGCTG')
    assert pair.offset == 59
    assert pair.head.name == 'seq1_901350_901788_1:0:0_0:0:0_21ca1/2'
    assert pair.swapped is False
    contig = kevlar.assemble.merge_pair(pair)
    print(contig)

    pair = kevlar.overlap.calc_offset(picorecord1, picorecord3,
                                      'TTTTTTGTTTCCCAAAGTAAGGCTG')
    assert pair.offset == 59
    assert pair.head.name == 'seq1_901350_901788_1:0:0_0:0:0_21ca1/2'
    assert pair.swapped is True
    newcontig = kevlar.assemble.merge_pair(pair)
    print(newcontig)
    assert kevlar.same_seq(contig, newcontig)
