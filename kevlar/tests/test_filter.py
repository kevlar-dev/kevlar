#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
import glob
import kevlar
from kevlar.count import load_sample_seqfile
from kevlar.seqio import AnnotatedReadSet as ReadSet
from kevlar.sequence import parse_augmented_fastx
from kevlar.tests import data_file
import khmer
import pytest
import sys
from tempfile import NamedTemporaryFile


def bogusrefr():
    maskfile = kevlar.tests.data_file('bogus-genome/refr.fa')
    return load_sample_seqfile([maskfile], 13, 1e7, count=False)


def test_alpha():
    readfile = data_file('collect.alpha.txt')
    filterer = kevlar.filter.filter(readfile, memory=500)
    validated = list(filterer)
    assert len(validated) == 8
    badkmers = ['CAGGCCAGGGATCGCCGTG']
    goodkmers = [
        'AGGGGCGTGACTTAATAAG', 'GGGCGTGACTTAATAAGGT',
        'TAGGGGCGTGACTTAATAA', 'GGGGCGTGACTTAATAAGG',
    ]
    for record in  validated:
        for kmer in record.annotations:
            seq = record.ikmerseq(kmer)
            assert seq not in badkmers and kevlar.revcom(seq) not in badkmers
            assert seq in goodkmers or kevlar.revcom(seq) in goodkmers


@pytest.mark.parametrize('mask,nkmers,nkmerinstances',[
    (None, 424, 5782),
    (bogusrefr(), 424, 5782),
    (kevlar.sketch.load(data_file('bogus-genome/mask.nt')), 13, 171)
])
def test_ctrl3(mask, nkmers, nkmerinstances):
    readfile = data_file('trio1/novel_3_1,2.txt')
    ikmers = defaultdict(int)
    for read in kevlar.filter.filter(readfile, memory=1e7, mask=mask):
        for ikmer in read.annotations:
            kmerseq = kevlar.revcommin(read.ikmerseq(ikmer))
            ikmers[kmerseq] += 1
    assert len(ikmers) == nkmers
    assert sum(ikmers.values()) == nkmerinstances


def test_filter_abundfilt():
    readfile = data_file('worm.augfasta')
    ikmers = defaultdict(int)
    filt = kevlar.filter.filter(readfile, memory=1000, casemin=5, ctrlmax=0)
    validated = list(filt)
    assert len(validated) == 5

    for read in validated:
        for ikmer in read.annotations:
            kmerseq = kevlar.revcommin(read.ikmerseq(ikmer))
            ikmers[kmerseq] += 1
    assert len(ikmers) == 1
    assert sum(ikmers.values()) == 5


def test_filter_main(capsys):
    arglist = [
        'filter',
        '--mask', data_file('bogus-genome/mask.nt'),
        '--memory', '10M', '--max-fpr', '0.001', '--case-min', '6',
        kevlar.tests.data_file('trio1/novel_3_1,2.txt'),
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.filter.main(args)

    out, err = capsys.readouterr()
    assert 'Processed 178 reads' in err
    assert 'Validated 18 reads' in err
