#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import sys
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import khmer
import kevlar


@pytest.mark.parametrize('filename,creads,ckmers,ukmers,kmerinst', [
    ('tests/data/trio1/novel_1_1,2.txt', 18, 684, 13, 171),
    ('tests/data/trio1/novel_2_1,2.txt', 53, 2014, 39, 540),
    ('tests/data/trio1/novel_3_1,2.txt', 178, 6764, 424, 5782),
    ('tests/data/trio1/novel_4_1,2.txt', 16, 608, 11, 99),
    ('tests/data/trio1/novel_5_3,4.txt', 16, 553, 9, 109),
    ('tests/data/trio1/novel_6_5,6.txt', 9, 336, 2, 9),
])
def test_load_input(filename, creads, ckmers, ukmers, kmerinst):
    ng = khmer.Nodegraph(13, 3e5 / 4, 4)
    vs = kevlar.VariantSet()
    with open(filename, 'r') as infile:
        nreads, nkmers = kevlar.collect.load_input(infile, ng, vs)
    assert nreads == creads
    assert nkmers == ckmers
    assert vs.nkmers == ukmers
    assert vs.nkmerinst == kmerinst
