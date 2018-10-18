#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import kevlar
from kevlar.readgraph import ReadGraph
from kevlar.tests import data_file


@pytest.mark.parametrize('partfile,edges,strictedges', [
    ('connectivity-1311.augfastq', 30, 11),
    ('connectivity-1541.augfastq', 31, 12),
])
def test_populate(partfile, edges, strictedges):
    with kevlar.open(data_file(partfile), 'r') as fh:
        reader = kevlar.parse_augmented_fastx(fh)
        reads = list(reader)
    rg = ReadGraph()
    rg.load(reads)
    rg.populate_edges()
    assert rg.number_of_edges() == pytest.approx(edges, 1)
    rg = ReadGraph()
    rg.load(reads)
    rg.populate_edges(strict=True)
    assert rg.number_of_edges() == pytest.approx(strictedges, 1)
