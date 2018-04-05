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


def test_gap_center_aligned():
    query = next(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file('cigar/b.contig.fa'), 'r')
        )
    )
    target = next(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file('cigar/b.gdna.fa'), 'r')
        )
    )
    cigar, score = kevlar.align(target.sequence, query.sequence)
    tok = AlignmentTokenizer(query.sequence, target.sequence, cigar)
    assert len(tok.blocks) == 3
    assert tok._cigar == '41D150M50D'
    assert tok._origcigar == '41D144M50D6M'
