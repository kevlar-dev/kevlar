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
import khmer
import kevlar


@pytest.fixture
def bogusrefr():
    refr = khmer.Nodetable(13, 1e7 / 4, 4)
    refrfile = kevlar.tests.data_file('bogus-genome/refr.fa')
    refr.consume_seqfile(refrfile)
    return refr


@pytest.fixture
def contaminants():
    contam = khmer.Nodetable(13, 1e7 / 4, 4)
    contamfile = kevlar.tests.data_file('bogus-genome/contam1.fa')
    contam.consume_seqfile(contamfile)
    return contam


@pytest.fixture
def ctrl3():
    augfastq = kevlar.tests.data_file('trio1/novel_3_1,2.txt')
    readset, countgraph = kevlar.filter.load_input([augfastq], 13, 1e7)
    return readset, countgraph


def test_load_refr():
    infile = kevlar.tests.data_file('bogus-genome/refr.fa')
    refr = kevlar.filter.load_refr(infile, 25, 1e7)
    assert refr.get('GGCCCCGAACTAGGGGGCCTACGTT') > 0
    assert refr.get('GCTGGCTAAATTTTCATACTAACTA') > 0
    assert refr.get('G' * 25) == 0


def test_load_input():
    filelist = kevlar.tests.data_glob('collect.beta.?.txt')
    readset, countgraph = kevlar.filter.load_input(filelist, 19, 1e3)

    assert len(readset) == 8
    assert readset
    kmers = [
        'AGGGGCGTGACTTAATAAG', 'GGGCGTGACTTAATAAGGT',
        'TAGGGGCGTGACTTAATAA', 'GGGGCGTGACTTAATAAGG',
    ]
    for kmer in kmers:
        assert countgraph.get(kmer) == 8


def test_validate():
    filelist = kevlar.tests.data_glob('collect.alpha.txt')
    readset, countgraph = kevlar.filter.load_input(filelist, 19, 5e3)
    kevlar.filter.validate_and_print(readset, countgraph)

    assert readset.valid == (4, 32)
    assert len(readset) == 9
    assert readset.discarded == 1

    badkmers = ['CAGGCCAGGGATCGCCGTG']
    goodkmers = [
        'AGGGGCGTGACTTAATAAG', 'GGGCGTGACTTAATAAGGT',
        'TAGGGGCGTGACTTAATAA', 'GGGGCGTGACTTAATAAGG',
    ]
    for record in readset:
        for kmer in record.ikmers:
            assert kmer.sequence not in badkmers and \
                kevlar.revcom(kmer.sequence) not in badkmers
            assert kmer.sequence in goodkmers or \
                kevlar.revcom(kmer.sequence) in goodkmers


def test_validate_minabund():
    filelist = kevlar.tests.data_glob('collect.beta.?.txt')
    readset, countgraph = kevlar.filter.load_input(filelist, 19, 5e3)
    kevlar.filter.validate_and_print(readset, countgraph)
    assert readset.valid == (4, 32)

    readset, countgraph = kevlar.filter.load_input(filelist, 19, 5e3)
    kevlar.filter.validate_and_print(readset, countgraph, minabund=9)
    assert readset.valid == (0, 0)


def test_validate_withrefr():
    kmer = 'AGGGGCGTGACTTAATAAG'
    refr = khmer.Nodetable(19, 1e3, 2)
    refr.add(kmer)

    filelist = kevlar.tests.data_glob('collect.beta.?.txt')
    readset, countgraph = kevlar.filter.load_input(filelist, 19, 5e3)
    kevlar.filter.validate_and_print(readset, countgraph, refr)
    assert readset.valid == (3, 24)
    for record in readset:
        for ikmer in record.ikmers:
            assert ikmer.sequence != kmer
            assert kevlar.revcom(ikmer.sequence) != kmer


def test_validate_withcontam():
    contam = khmer.Nodetable(19, 1e3, 2)
    contam.consume('CCAGCTGCAGGCCAGGGATCGCCGTGGGCGGACGCCCATACCGCGATAGC')
    filelist = kevlar.tests.data_glob('collect.gamma.txt')

    # First, without second pass
    readset, countgraph = kevlar.filter.load_input(filelist, 19, 5e3)
    kevlar.filter.validate_and_print(readset, countgraph, contam=contam,
                                     minabund=8, skip2=True)
    assert readset.valid == (5, 35)
    assert readset.lowabund == (0, 0)
    assert readset.discarded == 0
    assert readset.contam == 9

    # Then, with second pass
    readset, countgraph = kevlar.filter.load_input(filelist, 19, 5e3)
    kevlar.filter.validate_and_print(readset, countgraph, contam=contam,
                                     minabund=8)
    assert readset.valid == (4, 32)
    assert readset.lowabund == (2, 21)
    assert readset.discarded == 12
    assert readset.contam == 9


def test_ctrl3(ctrl3):
    readset, countgraph = ctrl3
    kevlar.filter.validate_and_print(readset, countgraph, minabund=6)
    assert readset.valid == (424, 5782)


def test_ctrl3_refr(ctrl3, bogusrefr):
    readset, countgraph = ctrl3
    kevlar.filter.validate_and_print(readset, countgraph, refr=bogusrefr,
                                     minabund=6)
    assert readset.valid == (424, 5782)


def test_ctrl3_refr_contam(ctrl3, bogusrefr, contaminants):
    readset, countgraph = ctrl3
    kevlar.filter.validate_and_print(readset, countgraph, refr=bogusrefr,
                                     contam=contaminants, minabund=6)
    assert readset.valid == (13, 171)
