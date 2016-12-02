#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
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
import kevlar


@pytest.fixture
def bogusargs():
    args = type('', (), {})()
    args.seqid = None
    args.genomemask = None
    args.out = StringIO()
    args.logfile = sys.stderr
    args.refr = 'tests/data/bogus-genome-refr.fa'
    args.reads = 'tests/data/bogus-reads.bam'
    return args


def test_basic(bogusargs):
    kevlar.dump.main(bogusargs)
    outputlines = bogusargs.out.getvalue().strip().split('\n')
    assert len(outputlines) == 5 * 4  # 5 records, 4 lines per record
    assert 'read2' in outputlines[0]
    assert 'read4' in outputlines[4]
    assert 'read6' in outputlines[8]
    assert 'read7' in outputlines[12]
    assert 'read8' in outputlines[16]


def test_seqid_filter(bogusargs):
    bogusargs.seqid = 'bogus-genome-chr1'
    kevlar.dump.main(bogusargs)
    outputlines = bogusargs.out.getvalue().strip().split('\n')
    assert len(outputlines) == 1 * 4  # 1 record, 4 lines per record
    assert 'read2' in outputlines[0]


def test_genomemask_filter(bogusargs):
    bogusargs.genomemask = 'tests/data/bogus-genome-refr-without-chr1.fa'
    bogusargs.maskmemory = 5e7
    bogusargs.mask_k = 13
    kevlar.dump.main(bogusargs)
    outputlines = bogusargs.out.getvalue().strip().split('\n')
    assert len(outputlines) == 3 * 4  # 3 records, 4 lines per record
    assert 'read2' in outputlines[0]
    assert 'read7' in outputlines[4]
    assert 'read8' in outputlines[8]


def test_seqid_genomemask_filters(bogusargs):
    bogusargs.seqid = 'bogus-genome-chr1'
    bogusargs.genomemask = 'tests/data/bogus-genome-refr-without-chr1.fa'
    bogusargs.maskmemory = 5e7
    bogusargs.mask_k = 13
    kevlar.dump.main(bogusargs)
    outputlines = bogusargs.out.getvalue().strip().split('\n')
    assert len(outputlines) == 3 * 4  # 3 records, 4 lines per record
    assert 'read2' in outputlines[0]
    assert 'read7' in outputlines[4]
    assert 'read8' in outputlines[8]


def test_indels(bogusargs):
    bogusargs.genomemask = 'tests/data/bogus-genome-refr-without-chr2.fa'
    bogusargs.maskmemory = 5e7
    bogusargs.mask_k = 13
    bogusargs.reads = 'tests/data/bogus-reads-indels.bam'
    kevlar.dump.main(bogusargs)
    outputlines = bogusargs.out.getvalue().strip().split('\n')
    assert len(outputlines) == 2 * 4  # 2 records, 4 lines per record
    assert 'read2' in outputlines[0]
    assert 'read4' in outputlines[4]
