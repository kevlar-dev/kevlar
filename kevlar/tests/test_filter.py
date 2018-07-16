#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import glob
import pytest
import sys
from tempfile import NamedTemporaryFile
import khmer
import kevlar
from kevlar.seqio import AnnotatedReadSet as ReadSet
from kevlar.sequence import parse_augmented_fastx


@pytest.fixture
def bogusrefr():
    mask = khmer.Nodetable(13, 1e7 / 4, 4)
    maskfile = kevlar.tests.data_file('bogus-genome/refr.fa')
    return kevlar.filter.load_mask([maskfile], 13, 1e7)


@pytest.fixture
def bogusrefrcontam():
    maskfile = kevlar.tests.data_file('bogus-genome/mask.nt')
    return kevlar.filter.load_mask([maskfile], 1, 1)


def test_load_mask():
    infile = kevlar.tests.data_file('bogus-genome/refr.fa')
    mask = kevlar.filter.load_mask([infile], 25, 1e7)
    assert mask.get('GGCCCCGAACTAGGGGGCCTACGTT') > 0
    assert mask.get('GCTGGCTAAATTTTCATACTAACTA') > 0
    assert mask.get('G' * 25) == 0

    assert kevlar.filter.load_mask(None, 25, 1e7) is None


def test_load_mask_multi_file():
    infiles = [
        kevlar.tests.data_file('bogus-genome/refr.fa'),
        kevlar.tests.data_file('bogus-genome/contam1.fa')
    ]
    mask = kevlar.filter.load_mask(infiles, 25, 1e7)
    assert mask.get('GGCCCCGAACTAGGGGGCCTACGTT') > 0  # reference
    assert mask.get('GCTGGCTAAATTTTCATACTAACTA') > 0  # reference
    assert mask.get('AATGTAGGTAGTTTTGTGCACAGTT') > 0  # contam
    assert mask.get('TCGCGCGCGTCCAAGTCGAGACCGC') > 0  # contam
    assert mask.get('G' * 25) == 0


def test_load_mask_save():
    infiles = [
        kevlar.tests.data_file('bogus-genome/refr.fa'),
        kevlar.tests.data_file('bogus-genome/contam1.fa')
    ]

    with NamedTemporaryFile(suffix='.nt') as table:
        mask = kevlar.filter.load_mask(infiles, 25, 1e7, savefile=table.name)
        newmask = khmer.Nodetable.load(table.name)
    assert newmask.get('GGCCCCGAACTAGGGGGCCTACGTT') > 0  # reference
    assert newmask.get('GCTGGCTAAATTTTCATACTAACTA') > 0  # reference
    assert newmask.get('AATGTAGGTAGTTTTGTGCACAGTT') > 0  # contam
    assert newmask.get('TCGCGCGCGTCCAAGTCGAGACCGC') > 0  # contam
    assert newmask.get('G' * 25) == 0


def test_load_mask_too_small():
    infile = kevlar.tests.data_file('bogus-genome/refr.fa')
    with pytest.raises(SystemExit) as se:
        mask = kevlar.filter.load_mask([infile], 25, 1e3)
    assert 'FPR too high, bailing out' in str(se)


def test_load_readset():
    filelist = kevlar.tests.data_glob('collect.beta.?.txt')
    readset = ReadSet(19, 1e3)
    for record in kevlar.seqio.afxstream(filelist):
        readset.add(record)

    assert len(readset) == 8
    assert readset
    kmers = [
        'AGGGGCGTGACTTAATAAG', 'GGGCGTGACTTAATAAGGT',
        'TAGGGGCGTGACTTAATAA', 'GGGGCGTGACTTAATAAGG',
    ]
    for kmer in kmers:
        assert readset._counts.get(kmer) == 8


def test_validate():
    filelist = kevlar.tests.data_glob('collect.alpha.txt')
    readset = ReadSet(19, 5e3)
    for record in kevlar.seqio.afxstream(filelist):
        readset.add(record)
    readset.validate()

    assert readset.valid == (4, 32)
    assert len(readset) == 9
    assert readset.discarded == 1

    badkmers = ['CAGGCCAGGGATCGCCGTG']
    goodkmers = [
        'AGGGGCGTGACTTAATAAG', 'GGGCGTGACTTAATAAGGT',
        'TAGGGGCGTGACTTAATAA', 'GGGGCGTGACTTAATAAGG',
    ]
    for record in readset:
        for kmer in record.annotations:
            seq = record.ikmerseq(kmer)
            assert seq not in badkmers and kevlar.revcom(seq) not in badkmers
            assert seq in goodkmers or kevlar.revcom(seq) in goodkmers


def test_validate_minabund():
    filelist = kevlar.tests.data_glob('collect.beta.?.txt')
    readset = ReadSet(19, 5e3)
    for record in kevlar.seqio.afxstream(filelist):
        readset.add(record)
    readset.validate()
    assert readset.valid == (4, 32)

    readset = ReadSet(19, 5e3)
    for record in kevlar.seqio.afxstream(filelist):
        readset.add(record)
    readset.validate(casemin=9)
    assert readset.valid == (0, 0)


def test_validate_with_mask():
    kmer = 'AGGGGCGTGACTTAATAAG'
    mask = khmer.Nodetable(19, 1e3, 2)
    mask.add(kmer)

    filelist = kevlar.tests.data_glob('collect.beta.?.txt')
    readset = ReadSet(19, 5e3, mask=mask)
    for record in kevlar.seqio.afxstream(filelist):
        readset.add(record)
    readset.validate()
    assert readset.valid == (3, 24)
    for record in readset:
        for ikmer in record.annotations:
            seq = record.ikmerseq(ikmer)
            assert seq != kmer
            assert kevlar.revcom(seq) != kmer


def test_ctrl3():
    augfastq = kevlar.tests.data_file('trio1/novel_3_1,2.txt')
    readset = ReadSet(13, 1e7)
    for record in parse_augmented_fastx(kevlar.open(augfastq, 'r')):
        readset.add(record)
    readset.validate(casemin=6)
    assert readset.valid == (424, 5782)


def test_ctrl3_refr(bogusrefr):
    augfastq = kevlar.tests.data_file('trio1/novel_3_1,2.txt')
    readset = ReadSet(13, 1e7, mask=bogusrefr)
    for record in parse_augmented_fastx(kevlar.open(augfastq, 'r')):
        readset.add(record)
    readset.validate()
    assert readset.valid == (424, 5782)


def test_ctrl3_refr_contam(bogusrefrcontam):
    augfastq = kevlar.tests.data_file('trio1/novel_3_1,2.txt')
    readset = ReadSet(13, 1e7, mask=bogusrefrcontam)
    for record in parse_augmented_fastx(kevlar.open(augfastq, 'r')):
        readset.add(record)
    readset.validate()
    assert readset.valid == (13, 171)


def test_filter_abundfilt():
    readset = ReadSet(31, 1000)
    augfastq = kevlar.tests.data_file('worm.augfasta')
    for record in parse_augmented_fastx(kevlar.open(augfastq, 'r')):
        readset.add(record)
    readset.validate(casemin=5, ctrlmax=0)
    assert readset.valid == (1, 5)
    assert readset.discarded == 2


def test_filter_main(capsys):
    arglist = [
        'filter',
        '--mask',
        kevlar.tests.data_file('bogus-genome/refr.fa'),
        kevlar.tests.data_file('bogus-genome/contam1.fa'),
        '--mask-memory', '10M',
        '--mask-max-fpr', '0.001',
        '--abund-memory', '10M',
        '--abund-max-fpr', '0.001',
        '--case-min', '6',
        '--ksize', '13',
        kevlar.tests.data_file('trio1/novel_3_1,2.txt'),
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.filter.main(args)

    out, err = capsys.readouterr()
    assert '171 instances of 13 distinct k-mers validated as novel' in err
