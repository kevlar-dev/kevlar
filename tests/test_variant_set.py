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


@pytest.fixture
def basicvset():
    vs = VariantSet()
    vs.add_kmer('TTTTT', 'read1', 'TTAAAAATT')
    vs.add_kmer('TTTTT', 'read2', 'TTAAAAATT')
    vs.add_kmer('TTTTT', 'read3', 'TTAAAAATT')
    vs.add_kmer('TTTTA', 'read4', 'TTAAAAATT')
    return vs


def test_vset_basic(basicvset):
    assert sorted(list(basicvset.kmers.keys())) == ['AAAAA', 'TAAAA']
    assert sorted(list(basicvset.read_ids.keys())) == \
        ['read1', 'read2', 'read3', 'read4']
    assert list(basicvset.linear_paths.keys()) == ['AATTTTTAA']


def test_vset_tables(basicvset):
    kout = StringIO()
    basicvset.kmer_table(kout)
    kouttest = ('AAAAA,TTTTT\tAATTTTTAA\tread1,read2,read3\n'
                'TAAAA,TTTTA\tAATTTTTAA\tread4\n')
    assert kout.getvalue() == kouttest

    pout = StringIO()
    basicvset.path_table(pout)
    pouttest = 'AATTTTTAA,TTAAAAATT\tread1,read2,read3,read4\tAAAAA,TAAAA\n'
    assert pout.getvalue() == pouttest

    rout = StringIO()
    basicvset.read_table(rout)
    routtest = ('read1\tAATTTTTAA\tAAAAA\n'
                'read2\tAATTTTTAA\tAAAAA\n'
                'read3\tAATTTTTAA\tAAAAA\n'
                'read4\tAATTTTTAA\tTAAAA\n')
    assert rout.getvalue() == routtest
