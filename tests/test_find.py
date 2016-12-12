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
def trio():
    readout = StringIO()
    kmerout = StringIO()
    pathout = StringIO()
    serrout = StringIO()

    args = type('', (), {})()
    args.controls = glob.glob('tests/data/trio1/ctrl?.counts')
    args.ctrl_max = 0
    args.case_min = 8
    args.out = readout
    args.flush = False
    args.kmers_out = kmerout
    args.paths_out = pathout
    args.logfile = serrout
    args.upint = 1000
    args.case = None
    args.case_fastq = None
    return args, pathout


def test_find_case1(trio):
    args, pathout = trio
    args.case = 'tests/data/trio1/case1.counts'
    args.case_fastq = 'tests/data/trio1/case1.fq'
    kevlar.find.main(args)

    path = ('GGCTATGGCGGAAGGGCACACCTAACCGCACCATTTGCCGTGGAAGCATAAAGGTCATCATTGAG'
            'GTGGTTCGTTCCGATACAGA')
    pathdata = pathout.getvalue()
    assert len(pathdata.strip().split('\n')) == 1
    assert pathdata.startswith(path)


def test_find_case2(trio):
    args, pathout = trio
    args.case = 'tests/data/trio1/case2.counts'
    args.case_fastq = 'tests/data/trio1/case2.fq'
    kevlar.find.main(args)

    pathdata = pathout.getvalue()
    assert len(pathdata.strip().split('\n')) == 4

    paths = sorted([ln.split(',')[0] for ln in pathdata.strip().split('\n')])
    assert paths == [
        'AAGGTAGTTCTCGGGGACCCTTAACGCACTTTAACCTTGATGCAGGT',
        'ACCAGGGGAGGTGAGAGTCAACCTTAGAACCGACCCATCCGTACGTAGCGATAGC',
        'ACCTGCATCAAGGTTAAAGTGCGTTAAGGGTCCCCGAGAACTACCTTGCCTTGCC',
        'GGCTATGGCGGAAGGGCACACCTAACCGCACCATTTGCCGTGGAAGCATAAAGGTCATCATTGAGGTGG'
        'TTCGTTCCGATACAGA',
    ]


def test_find_case3(trio):
    args, pathout = trio
    args.case = 'tests/data/trio1/case3.counts'
    args.case_fastq = 'tests/data/trio1/case3.fq'
    kevlar.find.main(args)

    pathdata = pathout.getvalue()
    assert len(pathdata.strip().split('\n')) == 59
