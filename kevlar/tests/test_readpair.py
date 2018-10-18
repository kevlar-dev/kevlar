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
import kevlar
from kevlar.sequence import KmerOfInterest, Record
from kevlar.readpair import ReadPair


@pytest.fixture
def record1():
    return Record(
        name='read1',
        sequence='GCTGCACCGATGTACGCAAA',
        annotations=[KmerOfInterest(5, 14, [15, 0, 0])],
    )


@pytest.fixture
def record2():
    return Record(
        name='read2',
        sequence='ACGCAAAGCTATTTAAAACC',
        annotations=[
            KmerOfInterest(5, 1, [15, 0, 0]),
            KmerOfInterest(5, 14, [19, 1, 0]),
        ],
    )


@pytest.fixture
def record2a():
    return Record(
        name='read2',
        sequence='ACGCAAAGCTATTTACGCAA',
        annotations=[
            KmerOfInterest(5, 1, [15, 0, 0]),
            KmerOfInterest(5, 15, [15, 0, 0]),
        ],
    )


@pytest.fixture
def record3():
    # reverse complement of record2
    return Record(
        name='read3',
        sequence='GGTTTTAAATAGCTTTGCGT',
        annotations=[
            KmerOfInterest(5, 1, [19, 1, 0]),
            KmerOfInterest(5, 14, [15, 0, 0]),
        ],
    )


@pytest.fixture
def record4():
    # similar to record2 but with a single nucleotide mismatch
    return Record(
        name='read4',
        sequence='ACGCAATGCTATTTAAAACC',
        annotations=[
            KmerOfInterest(5, 1, [15, 0, 0]),
            KmerOfInterest(5, 14, [19, 1, 0]),
        ],
    )


@pytest.fixture
def record5():
    return Record(
        name='read5',
        sequence='CTCTTCCGGCAGTCACTGTCAAGAGAGGGTGAACT',
        annotations=[
            KmerOfInterest(7, 15, [12, 0, 0]),
            KmerOfInterest(7, 16, [13, 0, 0]),
        ],
    )


@pytest.fixture
def record6():
    return Record(
        name='read6',
        sequence='TCACTGTCAAGAGAGGCCTACGGATTCGGTTACTG',
        annotations=[
            KmerOfInterest(7, 3, [12, 0, 0]),
            KmerOfInterest(7, 4, [13, 0, 0]),
        ],
    )


@pytest.fixture
def picorecord1():
    return Record(
        name='seq1_901350_901788_1:0:0_0:0:0_21ca1/2',
        sequence=('GTTTTTTTTTTGTTTCCCAAAGTAAGGCTGAGTGAACAATATTTTCTCATAGTTTTGAC'
                  'AAAAACAAAGGAATCCTTAGTTATTAAACTCGGGAGTTTGA'),
        annotations=[
            KmerOfInterest(25, 5, [19, 0, 0]),
            KmerOfInterest(25, 6, [18, 1, 0]),
            KmerOfInterest(25, 7, [18, 1, 0]),
            KmerOfInterest(25, 8, [18, 0, 0]),
            KmerOfInterest(25, 9, [17, 0, 0]),
        ],
    )


@pytest.fixture
def picorecord2():
    return Record(
        name='seq1_901428_901847_3:0:0_0:0:0_87d/1',
        sequence=('TTACATTTATTCGTTTGTGCAGGCTGAGACCTCACTTCCAACTGTAATCCAAAAGCTTA'
                  'GTTTTTTTTTTGTTTCCCAAAGTAAGGCTGAGTGAACAATA'),
        annotations=[
            KmerOfInterest(25, 64, [19, 0, 0]),
            KmerOfInterest(25, 65, [18, 1, 0]),
            KmerOfInterest(25, 66, [18, 1, 0]),
            KmerOfInterest(25, 67, [18, 0, 0]),
            KmerOfInterest(25, 68, [17, 0, 0]),
        ],
    )


@pytest.fixture
def picorecord3():
    return Record(
        name='seq1_901428_901847_3:0:0_0:0:0_87d/1',
        sequence=('TATTGTTCACTCAGCCTTACTTTGGGAAACAAAAAAAAAACTAAGCTTTTGGATTACAG'
                  'TTGGAAGTGAGGTCTCAGCCTGCACAAACGAATAAATGTAA'),
        annotations=[
            KmerOfInterest(25, 11, [17, 0, 0]),
            KmerOfInterest(25, 10, [18, 0, 0]),
            KmerOfInterest(25, 9, [18, 1, 0]),
            KmerOfInterest(25, 8, [18, 1, 0]),
            KmerOfInterest(25, 7, [19, 0, 0]),
        ],
    )


@pytest.fixture
def picorecord4():
    return Record(
        name='seqname',
        sequence=('TGTTCACTCAGCCTTACTTTGGGAAACAAAAAAAAAACTAAGCTTTTGGATTACAGTTG'
                  'GAAGTGAGGTCTCAGCCTGCACAAACGAATAAATG'),
        annotations=[
            KmerOfInterest(25, 8, [17, 0, 0]),
            KmerOfInterest(25, 7, [18, 0, 0]),
            KmerOfInterest(25, 6, [18, 1, 0]),
            KmerOfInterest(25, 5, [18, 1, 0]),
            KmerOfInterest(25, 4, [19, 0, 0]),
        ],
    )


@pytest.mark.parametrize('read1,read2,sameorientation', [
    (record1(), record2(), True),
    (record2(), record1(), True),
    (record1(), record3(), False),
    (record3(), record1(), False),
])
def test_basic(read1, read2, sameorientation):
    """Pair of reads sharing an interesting k-mer.

    GCTGCACCGATGTACGCAAA
                  |||||
                 ACGCAAAGCTATTTAAAACC

    We should get the same answer regardless of the order in which the reads
    are given or their orientation.
    """
    pair = ReadPair(read1, read2, 'CGCAA')
    print(pair.head.name, pair.tail.name, file=sys.stderr)
    print(str(pair), file=sys.stderr)
    assert pair.overlap == 7
    assert pair.offset == 13
    assert pair.sameorient is sameorientation
    str_reprs = [
        ('GCTGCACCGATGTACGCAAA\n'
         '              |||||\n'
         '             ACGCAAAGCTATTTAAAACC'),
        ('GGTTTTAAATAGCTTTGCGT\n'
         '              |||||\n'
         '             TTTGCGTACATCGGTGCAGC'),
    ]
    assert str(pair) in str_reprs


def test_kmer_multi_copy(record1, record2a):
    """Read contains multiple copies of the same k-mer.

    Make sure kevlar doesn't try to compute offset of reads when one of them
    contains more than one copy of the interesting k-mer.

    GCTGCACCGATGTACGCAAA
                  |||||         |||||
                 ACGCAAAGCTATTTACGCAA
    """
    pair = ReadPair(record1, record2a, 'CGCAA')
    assert pair.incompatible


def test_mismatch(record1, record4):
    """Pair of reads sharing an interesting k-mer, but with mismatch.

    GCTGCACCGATGTACGCAAA
                  |||||X
                 ACGCAATGCTATTTAAAACC

    The interesting k-mer is simply a seed. Assembly requires that the
    entire overlap between the two reads matches exactly, so the mismatch here
    should return a null offset indicating that the pair of reads is
    incompatible despite sharing an interesting k-mer.
    """
    pair = ReadPair(record1, record4, 'CGCAA')
    assert pair.incompatible


def test_big_mismatch(record5, record6):
    """Pair of reads sharing an interesting k-mer, but with mismatches.

    CTCTTCCGGCAGTCACTGTCAAGAGAGGGTGAACT
                   |||||||
                    |||||||     XXXXXXX
                TCACTGTCAAGAGAGGCCTACGGATTCGGTTACTG

    Not just a single mismatch, but extensive differences making the reads
    incompatible for assembly.
    """
    for ikmer in ['CTGTCAA', 'TGTCAAG']:
        pair = ReadPair(record5, record6, ikmer)
        assert pair.incompatible


def test_pico(picorecord1, picorecord2, picorecord3):
    pair1 = ReadPair(picorecord1, picorecord2, 'TTTTTTGTTTCCCAAAGTAAGGCTG')
    assert pair1.offset == 59
    assert pair1.head.read.name == 'seq1_901428_901847_3:0:0_0:0:0_87d/1'
    print(pair1.mergedseq)

    pair2 = ReadPair(picorecord1, picorecord3, 'TTTTTTGTTTCCCAAAGTAAGGCTG')
    assert pair2.offset == 59
    assert pair2.head.read.name == 'seq1_901428_901847_3:0:0_0:0:0_87d/1'
    print(pair2.mergedseq)

    assert kevlar.same_seq(pair1.mergedseq, pair2.mergedseq)


def test_pico_contains(picorecord3, picorecord4):
    pair = ReadPair(picorecord3, picorecord4, 'CACTCAGCCTTACTTTGGGAAACAA')
    print(pair.mergedseq)
    assert kevlar.same_seq(pair.mergedseq, picorecord3.sequence)
