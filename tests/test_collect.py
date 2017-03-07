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


def test_load_mask():
    mask = kevlar.collect.load_mask('tests/data/bogus-genome/refr.fa', 25, 1e7)
    assert mask.get('GGCCCCGAACTAGGGGGCCTACGTT') > 0
    assert mask.get('GCTGGCTAAATTTTCATACTAACTA') > 0
    assert mask.get('G' * 25) == 0


def test_recalc_abund_beta():
    filelist = glob.glob('tests/data/collect.beta.?.txt')
    countgraph = kevlar.collect.recalc_abund(filelist, 19, 1e3)

    kmers = [
        'AGGGGCGTGACTTAATAAG', 'GGGCGTGACTTAATAAGGT',
        'TAGGGGCGTGACTTAATAA', 'GGGGCGTGACTTAATAAGG',
    ]
    for kmer in kmers:
        assert countgraph.get(kmer) == 8


def test_load_novel_kmers_alpha():
    filelist = ['tests/data/collect.alpha.txt']
    countgraph = kevlar.collect.recalc_abund(filelist, 19, 1e3)
    vs = kevlar.collect.load_novel_kmers(filelist, countgraph)
    assert vs.nkmers == 4
    assert vs.nreads == 8

    badkmers = ['CAGGCCAGGGATCGCCGTG']
    goodkmers = [
        'AGGGGCGTGACTTAATAAG', 'GGGCGTGACTTAATAAGGT',
        'TAGGGGCGTGACTTAATAA', 'GGGGCGTGACTTAATAAGG',
    ]
    for kmer in badkmers:
        assert kmer not in vs.kmers and kevlar.revcom(kmer) not in vs.kmers
    for kmer in goodkmers:
        assert kmer in vs.kmers or kevlar.revcom(kmer) in vs.kmers


def test_load_novel_kmers_beta():
    filelist = glob.glob('tests/data/collect.beta.?.txt')
    countgraph = kevlar.collect.recalc_abund(filelist, 19, 1e3)
    vs = kevlar.collect.load_novel_kmers(filelist, countgraph)
    assert vs.nkmers == 4
    assert vs.nreads == 8

    countgraph = kevlar.collect.recalc_abund(filelist, 19, 1e3)
    vs = kevlar.collect.load_novel_kmers(filelist, countgraph, minabund=9)
    assert vs.nkmers == 0
    assert vs.nreads == 0


def test_load_novel_kmers_withmask():
    kmer = 'AGGGGCGTGACTTAATAAG'
    filelist = glob.glob('tests/data/collect.beta.?.txt')
    countgraph = kevlar.collect.recalc_abund(filelist, 19, 1e3)
    mask = khmer.Nodetable(19, 1e3, 2)
    mask.add(kmer)
    vs = kevlar.collect.load_novel_kmers(filelist, countgraph, mask)
    assert vs.nkmers == 3
    assert kmer not in vs.kmers and kevlar.revcom(kmer) not in vs.kmers
