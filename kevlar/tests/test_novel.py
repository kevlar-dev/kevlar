#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import glob
import pytest
import re
from tempfile import NamedTemporaryFile
import screed
import kevlar
from khmer import Counttable


def test_novel_banding_args():
    with pytest.raises(ValueError) as ve:
        reads = list(kevlar.novel.novel(None, [], [], numbands=4))
    assert 'Must specify `numbands` and `band` together' in str(ve)

    with pytest.raises(ValueError) as ve:
        reads = list(kevlar.novel.novel(None, [], [], band=0))
    assert 'Must specify `numbands` and `band` together' in str(ve)

    with pytest.raises(ValueError) as ve:
        reads = list(kevlar.novel.novel(None, [], [], numbands=4, band=-1))
    assert '`band` must be a value between 0 and 3' in str(ve)


def test_cli():
    args = kevlar.cli.parser().parse_args([
        'novel', '--case', 'case1.fq', '--control', 'cntl1.fq', '--control',
        'cntl2.fq', '-k', '17',
    ])
    assert args.ksize == 17
    assert args.case_min == 5
    assert args.ctrl_max == 1
    assert args.num_bands is None
    assert args.band is None

    args = kevlar.cli.parser().parse_args([
        'novel', '--num-bands', '8', '--band', '1', '--case', 'case1.fq',
        '--control', 'cntl1.fq', '--control', 'cntl2.fq',
    ])
    assert args.ksize == 31
    assert args.case_min == 5
    assert args.ctrl_max == 1
    assert args.num_bands == 8
    assert args.band == 1

    with pytest.raises(ValueError) as ve:
        args = kevlar.cli.parser().parse_args([
            'novel', '--case', 'case1.fq', '--control', 'cntl1.fq',
            '--band', '1'
        ])
        kevlar.novel.main(args)
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


@pytest.mark.long
@pytest.mark.parametrize('case,ctrl,mem', [
    ('trio1/case1.fq', 'trio1/ctrl[1,2].fq', '500K'),
    ('trio1/case2.fq', 'trio1/ctrl[1,2].fq', '1M'),
    ('trio1/case3.fq', 'trio1/ctrl[1,2].fq', '1M'),
    ('trio1/case4.fq', 'trio1/ctrl[1,2].fq', '500K'),
    ('trio1/case5.fq', 'trio1/ctrl[3,4].fq', '1M'),
    ('trio1/case6.fq', 'trio1/ctrl[5,6].fq', '1M'),
])
def test_novel_single_mutation(case, ctrl, mem, capsys):
    from sys import stdout, stderr
    casestr = kevlar.tests.data_file(case)
    ctrls = kevlar.tests.data_glob(ctrl)
    arglist = ['novel', '--case', casestr, '--ksize', '13', '--case-min', '8',
               '--control', ctrls[0], '--control', ctrls[1],
               '--ctrl-max', '0', '--memory', mem]
    args = kevlar.cli.parser().parse_args(arglist)
    args.out = None
    args.err = stderr
    kevlar.novel.main(args)
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


def test_novel_two_cases(capsys):
    cases = kevlar.tests.data_glob('trio1/case6*.fq')
    controls = kevlar.tests.data_glob('trio1/ctrl[5,6].fq')
    with NamedTemporaryFile(suffix='.ct') as case1ct, \
            NamedTemporaryFile(suffix='.ct') as case2ct, \
            NamedTemporaryFile(suffix='.ct') as ctrl1ct, \
            NamedTemporaryFile(suffix='.ct') as ctrl2ct:
        arglist = ['count', '--ksize', '19', '--memory', '1e7',
                   '--case', case1ct.name, cases[0],
                   '--case', case2ct.name, cases[1],
                   '--control', ctrl1ct.name, controls[0],
                   '--control', ctrl2ct.name, controls[1]]
        print(arglist)
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.count.main(args)

        arglist = ['novel', '--ksize', '19', '--memory', '1e7',
                   '--ctrl-max', '1', '--case-min', '7',
                   '--case', cases[0], '--case', cases[1],
                   '--case-counts', case1ct.name, case2ct.name,
                   '--control-counts', ctrl1ct.name, ctrl2ct.name]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.novel.main(args)
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

    kevlar.print_augmented_fastx(record, stdout)
    out, err = capsys.readouterr()
    assert read in out


def test_iter_read_multi_file():
    infiles = kevlar.tests.data_glob('bogus-genome/mask-chr[1,2].fa')
    print(infiles)
    records = [r for r in kevlar.multi_file_iter_khmer(infiles)]
    assert len(records) == 4


def test_novel_abund_screen(capsys):
    case = kevlar.tests.data_file('screen-case.fa')
    ctrl = kevlar.tests.data_file('screen-ctrl.fa')
    arglist = ['novel', '--ksize', '25', '--ctrl-max', '1', '--case-min', '8',
               '--case', case, '--control', ctrl, '--abund-screen', '3']
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.novel.main(args)

    out, err = capsys.readouterr()
    assert '>seq_error' not in out


def test_skip_until(capsys):
    readname = 'bogus-genome-chr1_115_449_0:0:0_0:0:0_1f4/1'
    case = kevlar.tests.data_file('trio1/case1.fq')
    ctrls = kevlar.tests.data_glob('trio1/ctrl[1,2].fq')
    arglist = ['novel', '--ctrl-max', '0', '--case-min', '6',
               '--skip-until', readname, '--upint', '50',
               '--case', case, '--control', ctrls[0], '--control', ctrls[1]]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.novel.main(args)

    out, err = capsys.readouterr()
    message = ('Found read bogus-genome-chr1_115_449_0:0:0_0:0:0_1f4/1 '
               '(skipped 1001 reads)')
    assert message in err
    assert '29 unique novel kmers in 14 reads' in err
