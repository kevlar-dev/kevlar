#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
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
import kevlar


@pytest.fixture
def trio1():
    readout = StringIO()
    kmerout = StringIO()
    pathout = StringIO()
    serrout = StringIO()

    args = type('', (), {})()
    args.controls = glob.glob('tests/data/trio1-ctrl?.counts')
    args.case = 'tests/data/trio1-case1.counts'
    args.ctrl_max = 0
    args.case_min = 8
    args.out = readout
    args.kmers_out = kmerout
    args.paths_out = pathout
    args.logfile = serrout
    args.upint = 1000
    args.case_fastq = 'tests/data/trio1-case1.fq'

    kevlar.find.main(args)
    return readout.getvalue(), kmerout.getvalue(), pathout.getvalue()


def test_basic(trio1):
    path = 'GATGACCTTTATGCTTCCACGGCAAATGGTGCGGTTAGGTGTGCCCTTCCGCCATAGCC'
    pathdata = trio1[2]
    assert len(pathdata.strip().split('\n')) == 1
    assert pathdata.startswith(path)
