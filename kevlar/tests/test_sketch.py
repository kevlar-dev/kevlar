#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import khmer
import kevlar


@pytest.mark.parametrize('filename,testkmer', [
    ('test.countgraph', 'TGGAACCGGCAACGACGAAAA'),
    ('test.smallcountgraph', 'CTGTACTACAGCTACTACAGT'),
    ('test.counttable', 'CCTGATATCCGGAATCTTAGC'),
    ('test.smallcounttable', 'GGGCCCCCATCTCTATCTTGC'),
    ('test.nodegraph', 'GGGAACTTACCTGGGGGTGCG'),
    ('test.nodetable', 'CTGTTCGATATGAGGAATCTG'),
])
def test_sketch_load(filename, testkmer):
    infile = kevlar.tests.data_file(filename)
    sketch = kevlar.sketch.load(infile)
    assert sketch.get(testkmer) > 0
    assert sketch.get('GATTACA' * 3) == 0


def test_sketch_load_badfilename():
    infile = kevlar.tests.data_file('test.notasketchtype')
    with pytest.raises(kevlar.sketch.KevlarSketchTypeError) as kste:
        sketch = kevlar.sketch.load(infile)
    assert ('sketch type from filename ' + infile) in str(kste)


@pytest.mark.parametrize('count,smallcount', [
    (True, True),
    (True, False),
    (False, False),
])
def test_allocate_sketch_graphy(count, smallcount):
    sequence = 'AATCAACGCTTCTTAATAGGCATAGTGTCTCTGCTGCGCATGGACGTGCCATAGCCACTACT'
    kmer = 'GCATAGTGTCTCTGCTGCGCA'

    sketch = kevlar.sketch.allocate(21, 1e4, 4, count, True, smallcount)
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

    sketch = kevlar.sketch.allocate(35, 1e4, 4, count, False, smallcount)
    sketch.consume(sequence)
    sketch.get(kmer) == 1
    kmer_hash = sketch.hash(kmer)
    with pytest.raises(ValueError) as ve:
        _ = sketch.reverse_hash(kmer_hash)
    assert 'not implemented' in str(ve)


def test_autoload():
    infile = kevlar.tests.data_file('test.nodegraph')
    sketch1 = kevlar.sketch.autoload(infile)
    assert sketch1.get('GGGAACTTACCTGGGGGTGCG') > 0

    infile = kevlar.tests.data_file('simple-genome-case-reads.fa.gz')
    sketch2 = kevlar.sketch.autoload(infile, ksize=25, table_size=1e7)
    assert sketch2.get('AGCTCAGACACTGGCGGTCTCTCCT') > 0

    sketch3 = kevlar.sketch.autoload(infile, ksize=25, table_size=1e7,
                                     count=True, graph=True, num_bands=4,
                                     band=0)
    assert sketch3.get('CAGCTGACCCACCGACACATAGGTT') > 0
