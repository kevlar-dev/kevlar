#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import khmer
import kevlar


@pytest.mark.parametrize('filename,count,graph,smallcount,testkmer', [
    ('test.countgraph', True, True, False, 'TGGAACCGGCAACGACGAAAA'),
    ('test.smallcountgraph', True, True, True, 'CTGTACTACAGCTACTACAGT'),
    ('test.counttable', True, False, False, 'CCTGATATCCGGAATCTTAGC'),
    ('test.smallcounttable', True, False, True, 'GGGCCCCCATCTCTATCTTGC'),
    ('test.nodegraph', False, True, False, 'GGGAACTTACCTGGGGGTGCG'),
    ('test.nodetable', False, False, False, 'CTGTTCGATATGAGGAATCTG'),
])
def test_load_sketch(filename, count, graph, smallcount, testkmer):
    infile = kevlar.tests.data_file(filename)
    sketch = kevlar.load_sketch(infile, count, graph, smallcount)
    assert sketch.get(testkmer) > 0
    assert sketch.get('GATTACA' * 3) == 0


@pytest.mark.parametrize('count,smallcount', [
    (True, True),
    (True, False),
    (False, False),
])
def test_allocate_sketch_graphy(count, smallcount):
    sequence = 'AATCAACGCTTCTTAATAGGCATAGTGTCTCTGCTGCGCATGGACGTGCCATAGCCACTACT'
    kmer = 'GCATAGTGTCTCTGCTGCGCA'

    sketch = kevlar.allocate_sketch(21, 1e4, 4, count, True, smallcount)
    sketch.consume(sequence)
    sketch.get(kmer) == 1
    kmer_hash = sketch.hash(kmer)
    assert kevlar.same_seq(sketch.reverse_hash(kmer_hash), kmer)


@pytest.mark.parametrize('count,smallcount', [
    (True, True),
    (True, False),
    (False, False),
])
def test_allocate_sketch_non_graphy(count, smallcount):
    sequence = 'TGCCACGATCCGGCTATGGCGGAAGGGCACACCTAACCGCGATGACGGAGTAACTCGCAGCA'
    kmer = 'CTATGGCGGAAGGGCACACCTAACCGCGATGACGG'

    sketch = kevlar.allocate_sketch(35, 1e4, 4, count, False, smallcount)
    sketch.consume(sequence)
    sketch.get(kmer) == 1
    kmer_hash = sketch.hash(kmer)
    with pytest.raises(ValueError) as ve:
        _ = sketch.reverse_hash(kmer_hash)
    assert 'not implemented' in str(ve)
