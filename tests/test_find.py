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


def test_cli():
    args = kevlar.cli.parser().parse_args([
        'find', '--controls', 'cntl1.fq', 'cntl2.fq', '--cases', 'case1.fq',
        '-k', '17'
    ])
    assert args.ksize == 17
    assert args.case_min == 5
    assert args.ctrl_max == 1
    assert args.num_bands is None
    assert args.band is None

    args = kevlar.cli.parser().parse_args([
        'find', '--controls', 'cntl1.fq', 'cntl2.fq', '--cases', 'case1.fq',
        '--num-bands', '8', '--band', '1'
    ])
    assert args.ksize == 31
    assert args.case_min == 5
    assert args.ctrl_max == 1
    assert args.num_bands == 8
    assert args.band == 1

    with pytest.raises(ValueError) as ve:
        args = kevlar.cli.parser().parse_args([
            'find', '--controls', 'cntl1.fq', '--cases', 'case1.fq',
            '--band', '1'
        ])
        kevlar.find.main(args)
    assert 'Must specify --num-bands and --band together' in str(ve)


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
    ('case1', 'ctrl[1,2]', '500K'),
    ('case1', 'ctrl[1,2]', '1M'),
    ('case2', 'ctrl[1,2]', '1M'),
    ('case3', 'ctrl[1,2]', '1M'),
    ('case4', 'ctrl[1,2]', '500K'),
    ('case4', 'ctrl[1,2]', '1M'),
    ('case5', 'ctrl[3,4]', '1M'),
    ('case6', 'ctrl[5,6]', '1M'),
])
def test_find_single_mutation(case, ctrl, mem, capsys):
    from sys import stdout, stderr
    casestr = 'tests/data/trio1/{}.fq'.format(case)
    ctrls = glob.glob('tests/data/trio1/{}.fq'.format(ctrl))
    arglist = ['find', '--ksize', '13', '--case_min', '8', '--ctrl_max', '0',
               '--memory', mem, '--cases', casestr, '--controls'] + ctrls
    args = kevlar.cli.parser().parse_args(arglist)
    args.out = stdout
    args.err = stderr
    kevlar.find.main(args)
    out, err = capsys.readouterr()

    for line in out.split('\n'):
        if not line.endswith('#'):
            continue
        abundmatch = re.search('(\d+) (\d+) (\d+)#$', line)
        assert abundmatch, line
        case = int(abundmatch.group(1))
        ctl1 = int(abundmatch.group(2))
        ctl2 = int(abundmatch.group(3))
        assert case >= 8, line
        assert ctl1 == 0 and ctl2 == 0, line


def test_find_two_cases(capsys):
    from sys import stdout, stderr
    pattern = 'tests/data/trio1/case{}.fq'
    cases = [pattern.format(n) for n in ('6', '6b')]
    ctrls = glob.glob('tests/data/trio1/ctrl[5,6].fq')
    arglist = ['find', '--ksize', '19', '--memory', '1e7', '--ctrl_max', '1',
               '--case_min', '7']
    arglist += ['--cases'] + cases + ['--controls'] + ctrls
    args = kevlar.cli.parser().parse_args(arglist)
    args.out = stdout
    args.err = stderr
    kevlar.find.main(args)
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
