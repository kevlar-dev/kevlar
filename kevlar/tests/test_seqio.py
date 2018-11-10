#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from io import StringIO
import pytest
import kevlar
from kevlar.tests import data_file
from kevlar.seqio import AnnotatedReadSet as ReadSet
from kevlar.seqio import KevlarPartitionLabelError
from kevlar.sequence import KmerOfInterest, Record
import khmer
import pysam
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


def test_augfastx_reader():
    infilename = data_file('collect.beta.1.txt')
    infile = open(infilename, 'r')
    for n, record in enumerate(kevlar.parse_augmented_fastx(infile)):
        assert record.name.startswith('good')
        assert record.sequence == (
            'TTAACTCTAGATTAGGGGCGTGACTTAATAAGGTGTGGGCCTAAGCGTCT'
        )
        assert len(record.annotations) == 2
        for kmer in record.annotations:
            assert kmer.abund == (8, 0, 0)
    assert n == 7


def test_augfastx_reader_e1():
    infilename = data_file('example1.augfastq')
    infile = open(infilename, 'r')
    record = next(kevlar.parse_augmented_fastx(infile))

    assert record.name == 'e1'
    assert record.sequence == (
        'TTAACTCTAGATTAGGGGCGTGACTTAATAAGGTGTGGGCCTAAGCGTCT'
    )
    assert len(record.annotations) == 2

    ikmer = record.annotations[0]
    assert record.ikmerseq(ikmer) == 'AGGGGCGTGACTTAATAAG'
    assert ikmer.ksize == 19
    assert ikmer.offset == 13
    assert ikmer.abund == (12, 15, 1, 1)

    ikmer = record.annotations[1]
    assert record.ikmerseq(ikmer) == 'GGGCGTGACTTAATAAGGT'
    assert ikmer.ksize == 19
    assert ikmer.offset == 15
    assert ikmer.abund == (20, 28, 0, 1)


def test_augfastx_reader_e2():
    infilename = data_file('example2.augfastq')
    infile = open(infilename, 'r')
    record = next(kevlar.parse_augmented_fastx(infile))

    assert record.name == 'ERR894724.125497791/1'
    assert record.sequence == (
        'TAGCCAGTTTGGGTAATTTTAATTGTAAAACTTTTTTTTCTTTTTTTTTGATTTTTTTTTTTCAAGCAG'
        'AAGACGGCATACGAGCTCTTTTCACGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT'
    )
    assert len(record.annotations) == 2

    ikmer = record.annotations[0]
    assert record.ikmerseq(ikmer) == 'GGCATACGAGCTCTTTTCACGTGACTGGAGT'
    assert ikmer.ksize == 31
    assert ikmer.offset == 74
    assert ikmer.abund == (23, 0, 0)

    ikmer = record.annotations[1]
    assert record.ikmerseq(ikmer) == 'GCTCTTTTCACGTGACTGGAGTTCAGACGTG'
    assert ikmer.ksize == 31
    assert ikmer.offset == 83
    assert ikmer.abund == (23, 0, 0)


def test_augfastx_reader_withmates():
    instream = kevlar.open(data_file('seqs-mates.augfastq'), 'r')
    reader = kevlar.parse_augmented_fastx(instream)

    record = next(reader)
    assert len(record.annotations) == 5
    assert len(record.mates) == 1
    assert record.mates[0].startswith('CTGATAAGCAACTTCAGCAAA')

    record = next(reader)
    assert len(record.annotations) == 4
    assert len(record.mates) == 1
    assert record.mates[0].startswith('ATTAGAAAAAAAAAGTGCATT')

    record = next(reader)
    assert len(record.annotations) == 21
    assert len(record.mates) == 0

    record = next(reader)
    assert len(record.annotations) == 2
    assert record.mates[0].startswith('CAGATGTGTCTTGTGGGCAGT')

    with pytest.raises(StopIteration):
        next(reader)


def test_augfastx_writer():
    output = StringIO()
    record = Record(
        name='BasiliscusVulgarisRead84467/1',
        sequence='TTAACTCTAGATTAGGGGCGTGACTTAATAAGGTGTGGGCCTAAGCGTCT',
        quality='BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB',
        annotations=[
            KmerOfInterest(ksize=19, offset=13, abund=(12, 1, 1)),
            KmerOfInterest(ksize=19, offset=15, abund=(20, 0, 1)),
        ],
    )
    kevlar.print_augmented_fastx(record, output)
    record = Record(
        name='BasiliscusVulgarisRead90577/2',
        sequence='CTGTAATCCCAGCACTTTGGGAGGCCGAGGCAAGCAGATGATGCGGTCAG',
        quality='BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB',
        annotations=[
            KmerOfInterest(ksize=19, offset=1, abund=(5, 7, 9)),
            KmerOfInterest(ksize=19, offset=2, abund=(7, 10, 9)),
        ],
        mates=['CAGATGTGTCTTGTGGGCAGTGCAGCGGAGAGGTGCAAATATGGGTTTGG']
    )
    kevlar.print_augmented_fastx(record, output)
    record = Record(
        name='BasiliscusVulgarisRead99037/1',
        sequence='AGCACTTTGGGAGGCCGAGGCAAGCAGATGATGCGGTCAGGATTACAGAT',
        quality='BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB'
    )
    kevlar.print_augmented_fastx(record, output)

    assert output.getvalue() == """@BasiliscusVulgarisRead84467/1
TTAACTCTAGATTAGGGGCGTGACTTAATAAGGTGTGGGCCTAAGCGTCT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
             AGGGGCGTGACTTAATAAG          12 1 1#
               GGGCGTGACTTAATAAGGT          20 0 1#
@BasiliscusVulgarisRead90577/2
CTGTAATCCCAGCACTTTGGGAGGCCGAGGCAAGCAGATGATGCGGTCAG
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
 TGTAATCCCAGCACTTTGG          5 7 9#
  GTAATCCCAGCACTTTGGG          7 10 9#
#mateseq=CAGATGTGTCTTGTGGGCAGTGCAGCGGAGAGGTGCAAATATGGGTTTGG#
@BasiliscusVulgarisRead99037/1
AGCACTTTGGGAGGCCGAGGCAAGCAGATGATGCGGTCAGGATTACAGAT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
"""


def test_partition_reader_simple():
    infile = data_file('part-reads-simple.fa')
    readstream = kevlar.parse_augmented_fastx(kevlar.open(infile, 'r'))
    partitions = list(kevlar.parse_partitioned_reads(readstream))
    partitions = [part for partid, part in partitions]
    assert len(partitions) == 2
    assert len(partitions[0]) == 4
    assert len(partitions[1]) == 2


def test_partition_reader_mixed():
    infile = data_file('part-reads-mixed.fa')
    readstream = kevlar.parse_augmented_fastx(kevlar.open(infile, 'r'))
    with pytest.raises(KevlarPartitionLabelError) as ple:
        partitions = list(kevlar.parse_partitioned_reads(readstream))
    assert 'with and without partition labels' in str(ple)


def test_parse_single_partition():
    infile = data_file('part-reads-simple.fa')

    readstream = kevlar.parse_augmented_fastx(kevlar.open(infile, 'r'))
    partitions = list(kevlar.parse_single_partition(readstream, '1'))
    partitions = [part for partid, part in partitions]
    assert len(partitions) == 1
    assert len(partitions[0]) == 4

    readstream = kevlar.parse_augmented_fastx(kevlar.open(infile, 'r'))
    partitions = list(kevlar.parse_single_partition(readstream, '2'))
    partitions = [part for partid, part in partitions]
    assert len(partitions) == 1
    assert len(partitions[0]) == 2

    readstream = kevlar.parse_augmented_fastx(kevlar.open(infile, 'r'))
    partitions = list(kevlar.parse_single_partition(readstream, 'alFrED'))
    partitions = [part for partid, part in partitions]
    assert partitions == []


def test_parse_single_partition_nonpartitioned_reads():
    infile = data_file('dup.augfastq')
    readstream = kevlar.parse_augmented_fastx(kevlar.open(infile, 'r'))
    partitions = list(kevlar.parse_single_partition(readstream, '42'))
    assert partitions == []


@pytest.mark.parametrize('basename', [
    ('example2.augfastq'),
    ('example2.augfastq.gz'),
])
def test_kevlar_open(basename):
    infilename = data_file(basename)
    infile = kevlar.open(infilename, 'r')
    record = next(kevlar.parse_augmented_fastx(infile))

    assert record.name == 'ERR894724.125497791/1'
    assert record.sequence == (
        'TAGCCAGTTTGGGTAATTTTAATTGTAAAACTTTTTTTTCTTTTTTTTTGATTTTTTTTTTTCAAGCAG'
        'AAGACGGCATACGAGCTCTTTTCACGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT'
    )
    assert len(record.annotations) == 2


def test_ikmer_abund_after_recalc():
    """Ensure interesting k-mer abundances are correct after recalculation.

    The interesting k-mer has an advertised abundance of 28, but a true
    abundance (in `counts`) of 10. The readset "validate" function should check
    and correct this.
    """
    read = Record(
        name='read1',
        sequence='AAGCAGGGGTCTACATTGTCCTCGGGACTCGAGATTTCTTCGCTGT',
        annotations=[KmerOfInterest(17, 13, (28, 0, 0))],
    )
    rs = ReadSet(17, 4e5)
    rs.add(read)

    seq = 'TTCGTTCCCGAAGCAGGGGTCTACATTGTCCTCGGGACTCGAGATTTCTTCGCTGTTCCGTCCTTCA'
    for _ in range(9):
        rs._counts.consume(seq)

    assert read.annotations[0].abund[0] == 28

    rs.validate(casemin=8)
    assert rs.valid == (1, 1)
    assert read.annotations[0].abund[0] == 10


def test_bam_paired_reader_5pairs():
    bam = pysam.AlignmentFile(data_file('ssc-10reads.bam'))
    reader = kevlar.dump.bam_paired_reader(bam)
    for n, (read1, read2) in enumerate(reader, 1):
        assert read1.query_name == read2.query_name
    assert n == 5


def test_bam_paired_reader_mixed():
    bam = pysam.AlignmentFile(data_file('ssc-8reads-mixed.bam'))
    reader = kevlar.dump.bam_paired_reader(bam)
    pairs = list(reader)
    assert len(pairs) == 5

    assert pairs[0][0].query_name == pairs[0][1].query_name
    assert pairs[1][0] is not None and pairs[1][1] is None
    assert pairs[2][0].query_name == pairs[2][1].query_name
    assert pairs[3][0] is not None and pairs[3][1] is None
    assert pairs[4][0].query_name == pairs[4][1].query_name


def test_ikmer_out_of_bounds():
    fh = kevlar.open(data_file('out-of-bounds.augfastq.gz'), 'r')
    reader = kevlar.parse_augmented_fastx(fh)
    with pytest.raises(AssertionError) as ae:
        list(reader)
    assert "('TACGACAGAC', 'TACGACAGACA')" in str(ae)
