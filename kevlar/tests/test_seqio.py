#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import kevlar
from kevlar import KmerOfInterest
from kevlar.seqio import AnnotatedReadSet as ReadSet
import khmer
import screed
import shutil


@pytest.fixture
def bogusseqs():
    seq = '>seq1\nACGT\n>seq2 yo\nGATTACA\nGATTACA\n>seq3\tdescrip\nATGATGTGA'
    return seq.split('\n')


def test_parse_fasta(bogusseqs):
    seqs = {name: seq for name, seq in kevlar.seqio.parse_fasta(bogusseqs)}
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
    infilename = kevlar.tests.data_file('collect.beta.1.txt')
    infile = open(infilename, 'r')
    for n, record in enumerate(kevlar.parse_augmented_fastx(infile)):
        assert record.name.startswith('good')
        assert record.sequence == (
            'TTAACTCTAGATTAGGGGCGTGACTTAATAAGGTGTGGGCCTAAGCGTCT'
        )
        assert len(record.ikmers) == 2
        for kmer in record.ikmers:
            assert kmer.abund == [8, 0, 0]
    assert n == 7


def test_aug_fastq_reader_e1():
    infilename = kevlar.tests.data_file('example1.augfastq')
    infile = open(infilename, 'r')
    record = next(kevlar.parse_augmented_fastx(infile))

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
    infilename = kevlar.tests.data_file('example2.augfastq')
    infile = open(infilename, 'r')
    record = next(kevlar.parse_augmented_fastx(infile))

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


@pytest.mark.parametrize('basename', [
    ('example2.augfastq'),
    ('example2.augfastq.gz'),
])
def test_kevlar_open(basename):
    infilename = kevlar.tests.data_file(basename)
    infile = kevlar.open(infilename, 'r')
    record = next(kevlar.parse_augmented_fastx(infile))

    assert record.name == 'ERR894724.125497791/1'
    assert record.sequence == (
        'TAGCCAGTTTGGGTAATTTTAATTGTAAAACTTTTTTTTCTTTTTTTTTGATTTTTTTTTTTCAAGCAG'
        'AAGACGGCATACGAGCTCTTTTCACGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT'
    )
    assert len(record.ikmers) == 2


def test_ikmer_abund_after_recalc():
    """
    Ensure interesting k-mer abundances are correct after recalculation.

    The interesting k-mer has an advertised abundance of 28, but a true
    abundance (in `counts`) of 10. The readset "validate" function should check
    and correct this.
    """
    read = screed.Record(
        name='read1',
        sequence='AAGCAGGGGTCTACATTGTCCTCGGGACTCGAGATTTCTTCGCTGT',
        ikmers=[KmerOfInterest('CATTGTCCTCGGGACTC', 13, [28, 0, 0])],
    )
    rs = ReadSet(17, 4e5)
    rs.add(read)

    seq = 'TTCGTTCCCGAAGCAGGGGTCTACATTGTCCTCGGGACTCGAGATTTCTTCGCTGTTCCGTCCTTCA'
    for _ in range(9):
        rs._counts.consume(seq)

    assert read.ikmers[0].abund[0] == 28

    rs.validate(minabund=8)
    assert rs.valid == (1, 1)
    assert read.ikmers[0].abund[0] == 10
