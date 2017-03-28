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
    for n, record in enumerate(kevlar.parse_augmented_fastq(infile)):
        assert record.name.startswith('good')
        assert record.sequence == (
            'TTAACTCTAGATTAGGGGCGTGACTTAATAAGGTGTGGGCCTAAGCGTCT'
        )
        assert len(record.ikmers) == 2
        for kmer in record.ikmers:
            assert kmer.abund == [8, 0, 0]
    assert n == 7


def test_aug_fastq_reader_e1():
    infile = open('tests/data/example1.augfastq', 'r')
    record = next(kevlar.parse_augmented_fastq(infile))

    assert record.name == 'e1'
    assert record.sequence == (
        'TTAACTCTAGATTAGGGGCGTGACTTAATAAGGTGTGGGCCTAAGCGTCT'
    )
    assert len(record.ikmers) == 2

    assert record.ikmers[0].sequence == 'AGGGGCGTGACTTAATAAG'
    assert record.ikmers[0].offset == 13
    assert record.ikmers[0].abund == [12, 15, 1, 1]

    assert record.ikmers[1].sequence == 'GGGCGTGACTTAATAAGGT'
    assert record.ikmers[1].offset == 15
    assert record.ikmers[1].abund == [20, 28, 0, 1]


def test_aug_fastq_reader_e2():
    infile = open('tests/data/example2.augfastq', 'r')
    record = next(kevlar.parse_augmented_fastq(infile))

    assert record.name == 'ERR894724.125497791/1'
    assert record.sequence == (
        'TAGCCAGTTTGGGTAATTTTAATTGTAAAACTTTTTTTTCTTTTTTTTTGATTTTTTTTTTTCAAGCAG'
        'AAGACGGCATACGAGCTCTTTTCACGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT'
    )
    assert len(record.ikmers) == 2

    assert record.ikmers[0].sequence == 'GGCATACGAGCTCTTTTCACGTGACTGGAGT'
    assert record.ikmers[0].offset == 74
    assert record.ikmers[0].abund == [23, 0, 0]

    assert record.ikmers[1].sequence == 'GCTCTTTTCACGTGACTGGAGTTCAGACGTG'
    assert record.ikmers[1].offset == 83
    assert record.ikmers[1].abund == [23, 0, 0]
