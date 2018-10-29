#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import kevlar
from kevlar.cigar import AlignmentBlock, AlignmentTokenizer
from kevlar.tests import data_file
import pytest


@pytest.mark.parametrize('contig,gdna', [
    (data_file('cigar/a.contig.fa'), data_file('cigar/a.gdna.fa')),
    (data_file('cigar/b.contig.fa'), data_file('cigar/b.gdna.fa')),
    (data_file('cigar/c.contig.fa'), data_file('cigar/c.gdna.fa')),
    (data_file('cigar/d.contig.fa'), data_file('cigar/d.gdna.fa')),
    (data_file('phony-snv-01.contig.fa'), data_file('phony-snv-01.gdna.fa')),
    (data_file('phony-snv-02.contig.fa'), data_file('phony-snv-02.gdna.fa')),
])
def test_blocks(contig, gdna):
    query = next(kevlar.parse_augmented_fastx(kevlar.open(contig, 'r')))
    target = next(kevlar.parse_augmented_fastx(kevlar.open(gdna, 'r')))
    cigar, score = kevlar.align(target.sequence, query.sequence)
    tok = AlignmentTokenizer(query.sequence, target.sequence, cigar)
    for block in tok.blocks:
        assert block.type in ('M', 'D', 'I')
        if block.type in ('M', 'D'):
            assert len(block.target) == block.length
        else:
            assert block.target is None
        if block.type in ('M', 'I'):
            assert len(block.query) == block.length
        else:
            assert block.query is None


def test_nomargin():
    qfile = kevlar.open(data_file('nomargin-r-indel-contigs.augfasta'), 'r')
    tfile = kevlar.open(data_file('nomargin-r-gdna.fa'), 'r')
    query = next(kevlar.parse_augmented_fastx(qfile))
    target = next(kevlar.parse_augmented_fastx(tfile))
    cigar, score = kevlar.align(target.sequence, query.sequence)
    tok = AlignmentTokenizer(query.sequence, target.sequence, cigar)
    assert tok._cigar == tok._origcigar


@pytest.mark.parametrize('contig,gdna,newcigar,origcigar,nblocks', [
    ('b.contig.fa', 'b.gdna.fa', '41D150M50D', '41D144M50D6M', 3),
    ('d.contig.fa', 'd.gdna.fa', '39D129M4D43M6D', '39D129M4D29M6D14M', 5),
])
def test_gap_center_aligned(contig, gdna, newcigar, origcigar, nblocks):
    qfile = kevlar.open(data_file('cigar/' + contig), 'r')
    tfile = kevlar.open(data_file('cigar/' + gdna), 'r')
    query = next(kevlar.parse_augmented_fastx(qfile))
    target = next(kevlar.parse_augmented_fastx(tfile))
    cigar, score = kevlar.align(target.sequence, query.sequence)
    tok = AlignmentTokenizer(query.sequence, target.sequence, cigar)
    assert len(tok.blocks) == nblocks
    assert tok._cigar == newcigar
    assert tok._origcigar == origcigar
