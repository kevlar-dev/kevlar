#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import kevlar
import re


@pytest.fixture
def bogusseqs():
    seq = '>seq1\nACGT\n>seq2 yo\nGATTACA\nGATTACA\n>seq3\tdescrip\nATGATGTGA'
    return seq.split('\n')


def test_parse(bogusseqs):
    seqs = {defline: seq for defline, seq in kevlar.fasta.parse(bogusseqs)}
    assert seqs == {
        '>seq1': 'ACGT',
        '>seq2 yo': 'GATTACAGATTACA',
        '>seq3\tdescrip': 'ATGATGTGA',
    }


def test_seq_dict(bogusseqs):
    d = kevlar.fasta.parse_seq_dict(bogusseqs)
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
