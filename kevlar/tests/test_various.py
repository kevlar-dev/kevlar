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
    ('test.smallcounttable', True, False, True, 'TTGGTATTAAGTAGCACTCGG'),
    ('test.nodegraph', False, True, False, 'GGGAACTTACCTGGGGGTGCG'),
    ('test.nodetable', False, False, False, 'GCTTACGCGCAAGCGCCCCCA'),
])
def test_load_sketch(filename, count, graph, smallcount, testkmer):
    infile = kevlar.tests.data_file(filename)
    sketch = kevlar.load_sketch(infile, count, graph, smallcount)
    assert sketch.get(testkmer) > 0
    assert sketch.get('GATTACA' * 3) == 0
