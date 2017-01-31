#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import glob
import pytest
import re
import sys
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import kevlar
from khmer import Counttable


@pytest.fixture
def trio_args():
    args = type('', (), {})()

    args.controls = glob.glob('tests/data/trio1/ctrl[1,2].fq')
    args.ctrl_max = 0
    args.case_min = 8
    args.ksize = 13
    args.memory = 1e6
    args.max_fpr = 0.2
    args.out = StringIO()
    args.flush = False
    args.collapse = False
    args.batch = None
    args.upint = 1000
    args.logfile = StringIO()
    args.cases = ['tests/data/trio1/case1.fq']

    return args


@pytest.mark.parametrize('kmer', [
    ('ACCGTACAA' * 3),
    ('TTATAATAG' * 3),
    ('CGAAAAATT' * 3),
])
def test_assumptions(kmer):
    ct = Counttable(27, 1e5, 2)
    kmer_rc = kevlar.revcom(kmer)
    assert ct.hash(kmer) == ct.hash(kmer_rc)
    assert ct.get_kmer_hashes(kmer)[0] == ct.get_kmer_hashes(kmer_rc)[0]


@pytest.mark.parametrize('case,ctrl,mem', [
    ('case1', 'ctrl[1,2]', 5e5),
    ('case1', 'ctrl[1,2]', 1e6),
    ('case2', 'ctrl[1,2]', 1e6),
    ('case3', 'ctrl[1,2]', 1e6),
    ('case4', 'ctrl[1,2]', 5e5),
    ('case4', 'ctrl[1,2]', 1e6),
    ('case5', 'ctrl[3,4]', 1e6),
    ('case6', 'ctrl[5,6]', 1e6),
])
def test_find_single_mutation(case, ctrl, mem, trio_args):
    trio_args.memory = mem
    trio_args.cases = ['tests/data/trio1/{}.fq'.format(case)]
    trio_args.controls = glob.glob('tests/data/trio1/{}.fq'.format(ctrl))
    kevlar.find.main(trio_args)

    for line in trio_args.out.getvalue().split('\n'):
        if not line.endswith('#'):
            continue
        abundmatch = re.search('(\d+) (\d+) (\d+)#$', line)
        assert abundmatch, line
        case = int(abundmatch.group(1))
        ctl1 = int(abundmatch.group(2))
        ctl2 = int(abundmatch.group(3))
        assert case >= 8
        assert ctl1 == 0 and ctl2 == 0


def test_kmer_rep_in_read():
    read = ('AGGATGAGGATGAGGATGAGGATGAGGATGAGGATGAGGATGAGGATGAGGATGAGGATGAGGAT'
            'GAGGATGAGGATGAGGAT')
    record = type('', (), {})()
    record.sequence = read
    record.name = 'reqseq'
    kmers = dict()
    kmers[2] = ['GATGAGGATGAGGATGAGGATGAGG', 11, 1, 0]
    kmers[8] = ['GATGAGGATGAGGATGAGGATGAGG', 11, 1, 0]
    outstream = StringIO()
    kevlar.find.print_interesting_read(record, kmers, outstream)
    assert read in outstream.getvalue()

    outstream = StringIO()
    kmers[14] = ['GATGAGGATGAGGATGAGGATGAGG', 11, 1, 0]
    kmers[20] = ['GATGAGGATGAGGATGAGGATGAGG', 11, 1, 0]
    kevlar.find.print_interesting_read(record, kmers, outstream)
    assert outstream.getvalue() == ''


def test_iter_screed():
    pattern = 'tests/data/bogus-genome/mask-chr{}.fa'
    infiles = [pattern.format(n) for n in (1, 2)]
    records = [r for r in kevlar.find.iter_screed(infiles)]
    assert len(records) == 4
