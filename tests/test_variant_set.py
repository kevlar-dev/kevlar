#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import kevlar
from kevlar.variantset import VariantSet


def test_vset_basic():
    vs = VariantSet()
    vs.add_kmer('TTTTT', 'read1', 'TTAAAAATT')
    vs.add_kmer('TTTTT', 'read2', 'TTAAAAATT')
    vs.add_kmer('TTTTT', 'read3', 'TTAAAAATT')
    vs.add_kmer('TTTTA', 'read4', 'TTAAAAATT')

    assert len(vs.kmers) == 2
    assert 'AAAAA' in vs.kmers
    assert 'TAAAA' in vs.kmers
