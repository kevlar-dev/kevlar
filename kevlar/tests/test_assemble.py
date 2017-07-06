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
from networkx import connected_components
import kevlar
from kevlar import KmerOfInterest
from kevlar.seqio import load_reads_and_kmers
from kevlar.assemble import (merge_pair, merge_and_reannotate)
from kevlar.overlap import (OverlappingReadPair, INCOMPATIBLE_PAIR,
                            calc_offset, graph_init_strict)
from kevlar.tests import data_file


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


@pytest.fixture
def record10():
    return screed.Record(
        name='read10',
        sequence=('CAGGTCCCCACCCGGATACTTGAAGCAGGCAGCCT'),
        ikmers=[
            KmerOfInterest('TCCCCACCCGGATACTT', 4, [28, 0, 0]),
            KmerOfInterest('CCCCACCCGGATACTTG', 5, [26, 0, 0]),
            KmerOfInterest('CCCGGATACTTGAAGCA', 10, [21, 0, 0]),
        ],
    )


@pytest.fixture
def record11():
    return screed.Record(
        name='read11',
        sequence='CCCGGATACTTGAAGCAGGCAGC',
        ikmers=[KmerOfInterest('CCCGGATACTTGAAGCA', 0, [21, 0, 0])],
    )


@pytest.fixture
def record12():
    return screed.Record(
        name='read12',
        sequence='CCCGGATACTTGAAGCAGGCAcC',
        ikmers=[KmerOfInterest('CCCGGATACTTGAAGCA', 0, [21, 0, 0])],
    )


def test_merge_pair(record1, record2, record4):
    """
    Assemble a compatible overlapping read pair.

    GCTGCACCGATGTACGCAAA
                  |||||                 -->   GCTGCACCGATGTACGCAAAGCTATTTAAAACC
                 ACGCAAAGCTATTTAAAACC
    """
    pair = OverlappingReadPair(tail=record1, head=record2, offset=13,
                               overlap=7, sameorient=True, swapped=False)
    assert merge_pair(pair) == 'GCTGCACCGATGTACGCAAAGCTATTTAAAACC'

    pair = OverlappingReadPair(tail=record1, head=record4, offset=13,
                               overlap=7, sameorient=True, swapped=False)
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
                               overlap=7, sameorient=True, swapped=False)
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
                               overlap=7, sameorient=False, swapped=False)
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
                               overlap=21, sameorient=True, swapped=False)
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
                               overlap=21, sameorient=False, swapped=False)
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


def test_merge_and_reannotate_contained(record7, record10):
    """
    Test merge/reannotation with containment.

    CAGGTCCCCACCCGGATACTTGAAGCAGGCAGCCTCAAGGTATGTGAGGCGATAACTCAA
        |||||||||||||||||
         |||||||||||||||||
              |||||||||||||||||
                                          |||||||||||||||||
                                           |||||||||||||||||
                                            |||||||||||||||||
    CAGGTCCCCACCCGGATACTTGAAGCAGGCAGCCT
        |||||||||||||||||
         |||||||||||||||||
              |||||||||||||||||
    """
    pair = OverlappingReadPair(tail=record7, head=record10, offset=0,
                               overlap=35, sameorient=True, swapped=False)
    newrecord = merge_and_reannotate(pair, 'ContainedAtOne')
    assert newrecord.name == 'ContainedAtOne'
    assert newrecord.sequence == record7.sequence
    assert newrecord.ikmers == record7.ikmers


def test_merge_and_reannotate_contained_with_offset(record7, record11):
    """
    Test merge/reannotation with containment and offset.

    CAGGTCCCCACCCGGATACTTGAAGCAGGCAGCCTCAAGGTATGTGAGGCGATAACTCAA
        |||||||||||||||||
         |||||||||||||||||
              |||||||||||||||||
                                          |||||||||||||||||
                                           |||||||||||||||||
                                            |||||||||||||||||
              CCCGGATACTTGAAGCAGGCAGC
              |||||||||||||||||
    """
    pair = OverlappingReadPair(tail=record7, head=record11, offset=10,
                               overlap=23, sameorient=True, swapped=False)
    newrecord = merge_and_reannotate(pair, 'ContainedAtOffset')
    assert newrecord.sequence == record7.sequence
    assert newrecord.ikmers == record7.ikmers


def test_merge_contained_with_offset_and_error(record7, record12):
    """
    Test merge with containment, offset, and a sequencing error.

    CAGGTCCCCACCCGGATACTTGAAGCAGGCAGCCTCAAGGTATGTGAGGCGATAACTCAA
        |||||||||||||||||
         |||||||||||||||||
              |||||||||||||||||
                                          |||||||||||||||||
                                           |||||||||||||||||
                                            |||||||||||||||||
              CCCGGATACTTGAAGCAGGCAcC
              |||||||||||||||||    *
    """
    pair = calc_offset(record7, record12, 'CCCGGATACTTGAAGCA')
    assert pair == INCOMPATIBLE_PAIR

    pair = OverlappingReadPair(tail=record7, head=record12, offset=10,
                               overlap=23, sameorient=True, swapped=False)
    with pytest.raises(AssertionError) as ae:
        newrecord = merge_and_reannotate(pair, 'WontWork')
    assert 'attempted to assemble incompatible reads' in str(ae)


def test_load_reads_and_kmers():
    """Make sure augmented records are loaded correctly."""
    instream = open(data_file('var1.reads.augfastq'), 'r')
    reads, kmers = load_reads_and_kmers(instream, logstream=None)
    assert len(reads) == 10
    assert len(kmers) == 7

    readname = 'read8f start=8,mutations=0'
    assert reads[readname].sequence == ('CACTGTCCTTACAGGTGGATAGTCGCTTTGTAATAAA'
                                        'AGAGTTACACCCCGGTTTTTAGAAGTCTCGACTTTAA'
                                        'GGAAGTGGGCCTACGGCGGAAGCCGT')

    testset = set([
        'read2f start=13,mutations=0', 'read8f start=8,mutations=0',
        'read10f start=34,mutations=0', 'read13f start=49,mutations=1',
        'read15f start=54,mutations=1', 'read16f start=13,mutations=1',
        'read22f start=5,mutations=0', 'read35f start=25,mutations=0',
        'read37f start=9,mutations=0'
    ])
    assert kmers['CCGGTTTTTAGAAGTCTCGACTTTAAGGA'] == testset


def test_graph_init():
    """Test graph initialization."""
    instream = open(data_file('var1.reads.augfastq'), 'r')
    reads, kmers = load_reads_and_kmers(instream, logstream=None)
    graph = graph_init_strict(reads, kmers, maxabund=None, logstream=None)

    # 10 reads in the file, but read16f has no valid connections due to error
    assert len(graph.nodes()) == 9

    # The given read shares its interesting k-mer and has compatible overlaps
    # with 6 other reads (read13f and read15f have errors).
    r23name = 'read23f start=67,mutations=0'
    assert len(graph[r23name]) == 6

    # Test the values of one of the edges.
    r35name = 'read35f start=25,mutations=0'
    assert graph[r23name][r35name]['offset'] == 42
    assert graph[r23name][r35name]['overlap'] == 58

    # Should all be a single CC
    assert len(list(connected_components(graph))) == 1

    r8name = 'read8f start=8,mutations=0'
    r37name = 'read37f start=9,mutations=0'
    assert graph[r37name][r8name]['offset'] == 1
    assert graph[r37name][r8name]['overlap'] == 99
    pair = OverlappingReadPair(
        tail=reads[r8name], head=reads[r37name], offset=1, overlap=99,
        sameorient=True, swapped=False
    )
    assert merge_pair(pair) == ('CACTGTCCTTACAGGTGGATAGTCGCTTTGTAATAAAAGAGTTAC'
                                'ACCCCGGTTTTTAGAAGTCTCGACTTTAAGGAAGTGGGCCTACGG'
                                'CGGAAGCCGTC')


def test_assembly_round2():
    instream = open(data_file('var1.round2.augfastq'), 'r')
    reads, kmers = load_reads_and_kmers(instream, logstream=None)
    contig, read = reads['contig1'], reads['read22f start=5,mutations=0']
    pair = calc_offset(contig, read, 'AAGTCTCGACTTTAAGGAAGTGGGCCTAC')
    assert pair.tail == read
    assert pair.head == contig
    assert merge_pair(pair) == ('TATCACTGTCCTTACAGGTGGATAGTCGCTTTGTAATAAAAGAGT'
                                'TACACCCCGGTTTTTAGAAGTCTCGACTTTAAGGAAGTGGGCCTA'
                                'CGGCGGAAGCCGTC')


def test_assembly_contigs():
    instream = open(data_file('AluContigs.augfastq'), 'r')
    reads, kmers = load_reads_and_kmers(instream, logstream=None)
    contig6, contig7 = reads['contig6'], reads['contig7']
    pair = calc_offset(contig6, contig7, 'AAAGTTTTCTTAAAAACATATATGGCCGGGC')
    assert pair.offset == 50
    assert pair.overlap == 85
    assert pair.tail == contig6
    newrecord = merge_and_reannotate(pair, 'newcontig')
    assert newrecord.sequence == ('TTGCCCAGGCTGGTCTCAAACTCCTGAGCTCAAAGCGATCTGT'
                                  'CGGCCTGGGCATCCAAAAAAAGTTTTCTTAAAAACATATATGG'
                                  'CCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAG'
                                  'GCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCC'
                                  'TGGCTAACACG')


def test_assemble_main(capsys):
    cliargs = ['assemble', data_file('var1.reads.augfastq')]
    args = kevlar.cli.parser().parse_args(cliargs)
    kevlar.assemble.main(args)
    out, err = capsys.readouterr()
    contig = ('TATCACTGTCCTTACAGGTGGATAGTCGCTTTGTAATAAAAGAGTTACACCCCGGTTTTTAGA'
              'AGTCTCGACTTTAAGGAAGTGGGCCTACGGCGGAAGCCGTCTCTAATGGACTCAAGGACCTGA'
              'ATCCAACTAGAGGAGCTTGCCAa')
    assert contig in out
