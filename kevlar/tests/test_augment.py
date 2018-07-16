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
import kevlar
from kevlar.augment import augment
from kevlar.tests import data_file


def test_augment_contigs():
    augfh = kevlar.open(data_file('snorkel.augfastq'), 'r')
    augreads = kevlar.parse_augmented_fastx(augfh)
    nakedfh = kevlar.open(data_file('snorkel-contig.fasta'), 'r')
    nakedseq = kevlar.parse_augmented_fastx(nakedfh)
    augseqs = list(augment(augreads, nakedseq))
    assert len(augseqs) == 1
    assert len(augseqs[0].annotations) == 3

    offsets = [k.offset for k in augseqs[0].annotations]
    assert offsets == [17, 20, 22]


def test_augment_reads(capsys):
    arglist = [
        'augment',
        data_file('reaugment.augfastq'),
        data_file('reaugment.fq')
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.augment.main(args)

    out, err = capsys.readouterr()
    testout = open(data_file('reaugment.out'), 'r').read()
    print(out)
    print(testout)
    assert out == testout


@pytest.mark.parametrize('docollapse,nummates', [
    (True, 26),
    (False, 0),
])
def test_augment_contig_mates(docollapse, nummates):
    augfh = kevlar.open(data_file('deadbeef.augfastq.gz'), 'r')
    augreads = kevlar.parse_augmented_fastx(augfh)
    nakedfh = kevlar.open(data_file('deadbeef.contig.fa'), 'r')
    nakedseq = kevlar.parse_augmented_fastx(nakedfh)
    contigs = list(augment(augreads, nakedseq, collapsemates=docollapse))
    assert len(contigs) == 1
    assert len(contigs[0].mates) == nummates
    assert len(contigs[0].annotations) == 74


def test_augment_reads_mates():
    augfh = kevlar.open(data_file('deadbeef.augfastq.gz'), 'r')
    augreads = list(kevlar.parse_augmented_fastx(augfh))
    nakedfh = kevlar.open(data_file('deadbeef.fq.gz'), 'r')
    nakedseq = kevlar.parse_augmented_fastx(nakedfh)
    newreads = list(augment(augreads, nakedseq, upint=5))
    for oldread, newread in zip(augreads, newreads):
        assert oldread.sequence == newread.sequence
        assert oldread.mates == newread.mates
        assert oldread.annotations == newread.annotations


def test_augment_cli(capsys):
    arglist = [
        'augment',
        data_file('snorkel.augfastq'),
        data_file('snorkel-contig.fasta')
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.augment.main(args)

    out, err = capsys.readouterr()
    assert out.strip() == """>contig1
AGGTCTTCGATGCTAGCATTTTTACGACAGACAAAAACAAGATTACATTCCAAAATACATACCGCGCC
                 ATTTTTACGAC          8 0 0#
                    TTTACGACAGA          11 0 0#
                      TACGACAGACA          9 0 0#"""
