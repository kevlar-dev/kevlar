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
from kevlar.tests import data_file


def test_augment_contigs():
    augfh = kevlar.open(data_file('snorkel.augfastq'), 'r')
    augreads = kevlar.parse_augmented_fastx(augfh)
    nakedseq = screed.open(data_file('snorkel-contig.fasta'))
    augmentor = kevlar.augment.augment(augreads, nakedseq)
    augseqs = list(augmentor)
    assert len(augseqs) == 1
    assert len(augseqs[0].ikmers) == 3

    offsets = [k.offset for k in augseqs[0].ikmers]
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
    assert out == testout


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
