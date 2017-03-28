#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import kevlar


@pytest.fixture
def bogusseqs():
    seq = '>seq1\nACGT\n>seq2 yo\nGATTACA\nGATTACA\n>seq3\tdescrip\nATGATGTGA'
    return seq.split('\n')


def test_parse_fasta(bogusseqs):
    seqs = {defline: seq for defline, seq in kevlar.parse_fasta(bogusseqs)}
    assert seqs == {
        '>seq1': 'ACGT',
        '>seq2 yo': 'GATTACAGATTACA',
        '>seq3\tdescrip': 'ATGATGTGA',
    }


def test_seq_dict(bogusseqs):
    d = kevlar.seqio.parse_seq_dict(bogusseqs)
    assert d == {
        'seq1': 'ACGT',
        'seq2': 'GATTACAGATTACA',
        'seq3': 'ATGATGTGA',
    }


def test_aug_fastq_reader():
    infile = open('tests/data/collect.beta.1.txt', 'r')
    for n, aug_record in enumerate(kevlar.parse_augmented_fastq(infile)):
        record, kmers = aug_record
        assert record.name.startswith('good')
        assert record.sequence == (
            'TTAACTCTAGATTAGGGGCGTGACTTAATAAGGTGTGGGCCTAAGCGTCT'
        )
        assert len(kmers) == 2
        for offset in kmers:
            kmer, abundances = kmers[offset]
            assert abundances == [8, 0, 0]
    assert n == 7


def test_aug_fastq_reader_e1():
    infile = open('tests/data/example1.augfastq', 'r')
    record, kmers = next(kevlar.parse_augmented_fastq(infile))

    assert record.name == 'e1'
    assert record.sequence == (
        'TTAACTCTAGATTAGGGGCGTGACTTAATAAGGTGTGGGCCTAAGCGTCT'
    )
    assert len(kmers) == 2

    assert 13 in kmers
    assert kmers[13][0] == 'AGGGGCGTGACTTAATAAG'
    assert kmers[13][1] == [12, 15, 1, 1]

    assert 15 in kmers
    assert kmers[15][0] == 'GGGCGTGACTTAATAAGGT'
    assert kmers[15][1] == [20, 28, 0, 1]


def test_aug_fastq_reader_e2():
    infile = open('tests/data/example2.augfastq', 'r')
    record, kmers = next(kevlar.parse_augmented_fastq(infile))

    assert record.name == 'ERR894724.125497791/1'
    assert record.sequence == (
        'TAGCCAGTTTGGGTAATTTTAATTGTAAAACTTTTTTTTCTTTTTTTTTGATTTTTTTTTTTCAAGCAG'
        'AAGACGGCATACGAGCTCTTTTCACGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT'
    )
    assert len(kmers) == 2

    assert 74 in kmers
    assert kmers[74][0] == 'GGCATACGAGCTCTTTTCACGTGACTGGAGT'
    assert kmers[74][1] == [23, 0, 0]

    assert 83 in kmers
    assert kmers[83][0] == 'GCTCTTTTCACGTGACTGGAGTTCAGACGTG'
    assert kmers[83][1] == [23, 0, 0]
