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
import screed
import kevlar
from khmer import Counttable


@pytest.fixture
def trio_args(capsys):
    from sys import stdout, stderr
    args = type('', (), {})()

    args.controls = glob.glob('tests/data/trio1/ctrl[1,2].fq')
    args.ctrl_max = 0
    args.case_min = 8
    args.ksize = 13
    args.memory = 1e6
    args.max_fpr = 0.2
    args.out = stdout
    args.flush = False
    args.collapse = False
    args.batch = None
    args.upint = 1000
    args.logfile = stderr
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
def test_find_single_mutation(case, ctrl, mem, trio_args, capsys):
    trio_args.memory = mem
    trio_args.cases = ['tests/data/trio1/{}.fq'.format(case)]
    trio_args.controls = glob.glob('tests/data/trio1/{}.fq'.format(ctrl))
    kevlar.find.main(trio_args)
    out, err = capsys.readouterr()

    for line in out.split('\n'):
        if not line.endswith('#'):
            continue
        abundmatch = re.search('(\d+) (\d+) (\d+)#$', line)
        assert abundmatch, line
        case = int(abundmatch.group(1))
        ctl1 = int(abundmatch.group(2))
        ctl2 = int(abundmatch.group(3))
        assert case >= 8
        assert ctl1 == 0 and ctl2 == 0


def test_find_two_cases(trio_args, capsys):
    from sys import stdout, stderr
    pattern = 'tests/data/trio1/case{}.fq'
    trio_args.cases = [pattern.format(n) for n in ('6', '6b')]
    trio_args.controls = glob.glob('tests/data/trio1/ctrl[5,6].fq')
    trio_args.ctrl_max = 1
    trio_args.case_min = 7
    trio_args.ksize = 19
    trio_args.memory = 1e7
    trio_args.out = stdout
    trio_args.err = stderr
    kevlar.find.main(trio_args)
    out, err = capsys.readouterr()

    assert out.strip() != ''
    for line in out.split('\n'):
        if not line.endswith('#'):
            continue
        abundmatch = re.search('(\d+) (\d+) (\d+) (\d+)#$', line)
        assert abundmatch, line
        case1 = int(abundmatch.group(1))
        case2 = int(abundmatch.group(2))
        ctl1 = int(abundmatch.group(3))
        ctl2 = int(abundmatch.group(4))
        assert case1 >= 7 and case2 >= 7
        assert ctl1 <= 1 and ctl2 <= 1


def test_kmer_rep_in_read(capsys):
    from sys import stdout
    read = ('AGGATGAGGATGAGGATGAGGATGAGGATGAGGATGAGGATGAGGATGAGGATGAGGATGAGGAT'
            'GAGGATGAGGATGAGGAT')
    record = screed.Record(name='reqseq', sequence=read, ikmers=list())

    k1 = kevlar.KmerOfInterest(
        sequence='GATGAGGATGAGGATGAGGATGAGG',
        offset=2,
        abund=[11, 1, 0]
    )
    k2 = kevlar.KmerOfInterest(
        sequence='GATGAGGATGAGGATGAGGATGAGG',
        offset=8,
        abund=[11, 1, 0]
    )
    record.ikmers.extend([k1, k2])

    kevlar.print_augmented_fastq(record, stdout)
    out, err = capsys.readouterr()
    assert read in out


def test_iter_screed():
    pattern = 'tests/data/bogus-genome/mask-chr{}.fa'
    infiles = [pattern.format(n) for n in (1, 2)]
    records = [r for r in kevlar.find.iter_screed(infiles)]
    assert len(records) == 4
