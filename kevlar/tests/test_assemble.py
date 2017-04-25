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
from kevlar.assemble import (merge_pair, merge_and_reannotate,
                             OverlappingReadPair)


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


@pytest.fixture
def record7():
    return screed.Record(
        name='read7',
        sequence=('CAGGTCCCCACCCGGATACTTGAAGCAGGCAGCCTCAAGGTATGTGAGGC'
                  'GATAACTCAA'),
        ikmers=[
            KmerOfInterest('TCCCCACCCGGATACTT', 4, [28, 0, 0]),
            KmerOfInterest('CCCCACCCGGATACTTG', 5, [26, 0, 0]),
            KmerOfInterest('CCCGGATACTTGAAGCA', 10, [21, 0, 0]),
            KmerOfInterest('GGTATGTGAGGCGATAA', 38, [14, 0, 0]),
            KmerOfInterest('GTATGTGAGGCGATAAC', 39, [15, 1, 0]),
            KmerOfInterest('TATGTGAGGCGATAACT', 40, [15, 1, 1]),
        ],
    )


@pytest.fixture
def record8():
    return screed.Record(
        name='read8',
        sequence=('GTATGTGAGGCGATAACTCAAGACCACGGGAGCTCACTTCGTTGACGCGA'
                  'GCGCCTTGCT'),
        ikmers=[
            KmerOfInterest('GTATGTGAGGCGATAAC', 0, [15, 1, 0]),
            KmerOfInterest('TATGTGAGGCGATAACT', 1, [15, 1, 1]),
            KmerOfInterest('TGACGCGAGCGCCTTGC', 42, [39, 0, 0]),
            KmerOfInterest('GACGCGAGCGCCTTGCT', 43, [25, 0, 0]),
        ],
    )


@pytest.fixture
def record9():
    return screed.Record(
        name='read9',
        sequence=('AGCAAGGCGCTCGCGTCAACGAAGTGAGCTCCCGTGGTCTTGAGTTATCG'
                  'CCTCACATAC'),
        ikmers=[
            KmerOfInterest('AGCAAGGCGCTCGCGTC', 0, [25, 0, 0]),
            KmerOfInterest('GCAAGGCGCTCGCGTCA', 1, [39, 0, 0]),
            KmerOfInterest('GTTATCGCCTCACATAC', 42, [15, 1, 1]),
            KmerOfInterest('AGTTATCGCCTCACATA', 43, [15, 1, 0]),
        ],
    )


def test_calc_offset_same_orientation(record1, record2):
    """
    Compute offset of reads sharing an interesting k-mer, same orientation.

    GCTGCACCGATGTACGCAAA
                  |||||
                 ACGCAAAGCTATTTAAAACC
    """
    pair = kevlar.assemble.calc_offset(record1, record2, 'CGCAA')
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
    pair = kevlar.assemble.calc_offset(record1, record3, 'CGCAA')
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
    pair = kevlar.assemble.calc_offset(record1, record4, 'CGCAA')
    assert pair == kevlar.assemble.IncompatiblePair


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
        pair = kevlar.assemble.calc_offset(record5, record6, ikmer)
        assert pair == kevlar.assemble.IncompatiblePair


def test_merge_pair(record1, record2, record4):
    """
    Assemble a compatible overlapping read pair.

    GCTGCACCGATGTACGCAAA
                  |||||                 -->   GCTGCACCGATGTACGCAAAGCTATTTAAAACC
                 ACGCAAAGCTATTTAAAACC
    """
    pair = OverlappingReadPair(tail=record1, head=record2, offset=13,
                               overlap=7, sameorient=True)
    assert merge_pair(pair) == 'GCTGCACCGATGTACGCAAAGCTATTTAAAACC'

    pair = OverlappingReadPair(tail=record1, head=record4, offset=13,
                               overlap=7, sameorient=True)
    with pytest.raises(AssertionError) as ae:
        contig = merge_pair(pair)
    assert 'attempted to assemble incompatible reads' in str(ae)


def test_merge_and_reannotate_same_orientation(record1, record2):
    """
    Assemble a read pair and re-annotate the associated interesting k-mers.

    GCTGCACCGATGTACGCAAA
                  |||||                 -->   GCTGCACCGATGTACGCAAAGCTATTTAAAACC
                 ACGCAAAGCTATTTAAAACC                       *****        *****
    """
    pair = OverlappingReadPair(tail=record1, head=record2, offset=13,
                               overlap=7, sameorient=True)
    newrecord = merge_and_reannotate(pair, 'contig1')
    assert newrecord.name == 'contig1'
    assert newrecord.sequence == 'GCTGCACCGATGTACGCAAAGCTATTTAAAACC'
    assert len(newrecord.ikmers) == 2
    assert newrecord.ikmers[0].offset == 14
    assert newrecord.ikmers[1].offset == 27


def test_merge_and_reannotate_opposite_orientation(record1, record3):
    """
    Assemble a read pair and re-annotate the associated interesting k-mers.

    GCTGCACCGATGTACGCAAA
                  |||||                 -->   GCTGCACCGATGTACGCAAAGCTATTTAAAACC
                 ACGCAAAGCTATTTAAAACC                       *****        *****
    """
    pair = OverlappingReadPair(tail=record1, head=record3, offset=13,
                               overlap=7, sameorient=False)
    newrecord = merge_and_reannotate(pair, 'contig1')
    assert newrecord.name == 'contig1'
    assert newrecord.sequence == 'GCTGCACCGATGTACGCAAAGCTATTTAAAACC'
    assert len(newrecord.ikmers) == 2
    assert newrecord.ikmers[0].offset == 14
    assert newrecord.ikmers[1].offset == 27


def test_merge_and_reannotate_edge_case_same_orientation(record7, record8):
    """
    Test merge/reannotation with edge cases to check for fencepost errors.

    CAGGTCCCCACCCGGATACTTGAAGCAGGCAGCCTCAAGGTATGTGAGGCGATAACTCAA
        |||||||||||||||||
         |||||||||||||||||
              |||||||||||||||||
                                          |||||||||||||||||
                                           |||||||||||||||||
                                            |||||||||||||||||
                                           GTATGTGAGGCGATAACTCAAGACCACGGGAGCTCACTTCGTTGACGCGAGCGCCTTGCT
                                           |||||||||||||||||
                                            |||||||||||||||||
                                                                                     |||||||||||||||||
                                                                                      |||||||||||||||||
    """  # noqa
    pair = OverlappingReadPair(tail=record7, head=record8, offset=39,
                               overlap=21, sameorient=True)
    newrecord = merge_and_reannotate(pair, 'SoMeCoNtIg')
    assert newrecord.name == 'SoMeCoNtIg'
    assert newrecord.sequence == ('CAGGTCCCCACCCGGATACTTGAAGCAGGCAGCCTCAAGGTAT'
                                  'GTGAGGCGATAACTCAAGACCACGGGAGCTCACTTCGTTGACG'
                                  'CGAGCGCCTTGCT')
    assert len(newrecord.ikmers) == 8

    testseqs = [
        'TCCCCACCCGGATACTT', 'CCCCACCCGGATACTTG', 'CCCGGATACTTGAAGCA',
        'GGTATGTGAGGCGATAA', 'GTATGTGAGGCGATAAC', 'TATGTGAGGCGATAACT',
    ]
    testoffsets = [4, 5, 10, 38, 39, 40, 81, 82]
    for kmer, seq, offset in zip(newrecord.ikmers, testseqs, testoffsets):
        assert kmer.sequence == seq
        assert kmer.offset == offset


def test_merge_and_reannotate_edge_case_opposite_orientation(record7, record9):
    """
    Test merge/reannotation with edge cases to check for fencepost errors.

    CAGGTCCCCACCCGGATACTTGAAGCAGGCAGCCTCAAGGTATGTGAGGCGATAACTCAA
        |||||||||||||||||
         |||||||||||||||||
              |||||||||||||||||
                                          |||||||||||||||||
                                           |||||||||||||||||
                                            |||||||||||||||||
               reverse complement  ----->  GTATGTGAGGCGATAACTCAAGACCACGGGAGCTCACTTCGTTGACGCGAGCGCCTTGCT
                                           |||||||||||||||||
                                            |||||||||||||||||
                                                                                     |||||||||||||||||
                                                                                      |||||||||||||||||
    """  # noqa
    pair = OverlappingReadPair(tail=record7, head=record9, offset=39,
                               overlap=21, sameorient=False)
    newrecord = merge_and_reannotate(pair, 'SoMeCoNtIg')
    assert newrecord.name == 'SoMeCoNtIg'
    assert newrecord.sequence == ('CAGGTCCCCACCCGGATACTTGAAGCAGGCAGCCTCAAGGTAT'
                                  'GTGAGGCGATAACTCAAGACCACGGGAGCTCACTTCGTTGACG'
                                  'CGAGCGCCTTGCT')
    assert len(newrecord.ikmers) == 8

    testseqs = [
        'TCCCCACCCGGATACTT', 'CCCCACCCGGATACTTG', 'CCCGGATACTTGAAGCA',
        'GGTATGTGAGGCGATAA', 'GTATGTGAGGCGATAAC', 'TATGTGAGGCGATAACT',
    ]
    testoffsets = [4, 5, 10, 38, 39, 40, 81, 82]
    for kmer, seq, offset in zip(newrecord.ikmers, testseqs, testoffsets):
        assert kmer.sequence == seq
        assert kmer.offset == offset
