#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import glob
import pytest
import sys
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import khmer
import kevlar


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
