#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import filecmp
import glob
import json
import pytest
import re
from tempfile import NamedTemporaryFile, mkdtemp
import screed
from shutil import rmtree
import kevlar
from kevlar.tests import data_file, data_glob
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
    assert args.case_min == 6
    assert args.ctrl_max == 1
    assert args.num_bands is None
    assert args.band is None

    args = kevlar.cli.parser().parse_args([
        'novel', '--num-bands', '8', '--band', '1', '--case', 'case1.fq',
        '--control', 'cntl1.fq', '--control', 'cntl2.fq',
    ])
    assert args.ksize == 31
    assert args.case_min == 6
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


@pytest.mark.parametrize('case,ctrl', [
    ('microtrios/trio-li-proband.fq.gz', 'microtrios/trio-li-??ther.fq.gz'),
    ('microtrios/trio-na-proband.fq.gz', 'microtrios/trio-na-??ther.fq.gz'),
    ('microtrios/trio-k-proband.fq.gz',  'microtrios/trio-k-??ther.fq.gz'),
])
def test_novel_single_mutation(case, ctrl, capsys):
    casestr = data_file(case)
    ctrls = kevlar.tests.data_glob(ctrl)
    arglist = ['novel', '--case', casestr, '--ksize', '25', '--case-min', '7',
               '--control', ctrls[0], '--control', ctrls[1],
               '--num-bands', '2', '--band', '2',
               '--ctrl-max', '0', '--memory', '500K']
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.novel.main(args)
    out, err = capsys.readouterr()

    for line in out.split('\n'):
        if not line.endswith('#') or line.startswith('#mateseq'):
            continue
        abundmatch = re.search(r'(\d+) (\d+) (\d+)#$', line)
        assert abundmatch, line
        case = int(abundmatch.group(1))
        ctl1 = int(abundmatch.group(2))
        ctl2 = int(abundmatch.group(3))
        assert case >= 7, line
        assert ctl1 == 0 and ctl2 == 0, line


def test_novel_two_cases(capsys):
    cases = kevlar.tests.data_glob('trio1/case6*.fq')
    controls = kevlar.tests.data_glob('trio1/ctrl[5,6].fq')
    with NamedTemporaryFile(suffix='.ct') as case1ct, \
            NamedTemporaryFile(suffix='.ct') as case2ct, \
            NamedTemporaryFile(suffix='.ct') as ctrl1ct, \
            NamedTemporaryFile(suffix='.ct') as ctrl2ct:
        counttables = [case1ct, case2ct, ctrl1ct, ctrl2ct]
        seqfiles = cases + controls
        for ct, seqfile in zip(counttables, seqfiles):
            arglist = ['count', '--ksize', '19', '--memory', '1e7', ct.name,
                       seqfile]
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
        if not line.endswith('#') or line.startswith('#mateseq'):
            continue
        abundmatch = re.search(r'(\d+) (\d+) (\d+) (\d+)#$', line)
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
    record = kevlar.sequence.Record(name='reqseq', sequence=read)
    record.annotate('GATGAGGATGAGGATGAGGATGAGG', 2, (11, 1, 0))
    record.annotate('GATGAGGATGAGGATGAGGATGAGG', 8, (11, 1, 0))

    kevlar.print_augmented_fastx(record, stdout)
    out, err = capsys.readouterr()
    assert read in out


def test_iter_read_multi_file():
    infiles = kevlar.tests.data_glob('bogus-genome/mask-chr[1,2].fa')
    print(infiles)
    records = [r for r in kevlar.multi_file_iter_khmer(infiles)]
    assert len(records) == 4


def test_novel_abund_screen(capsys):
    case = data_file('screen-case.fa')
    ctrl = data_file('screen-ctrl.fa')
    arglist = ['novel', '--ksize', '25', '--ctrl-max', '1', '--case-min', '8',
               '--case', case, '--control', ctrl, '--abund-screen', '3']
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.novel.main(args)

    out, err = capsys.readouterr()
    assert '>seq_error' not in out


def test_skip_until(capsys):
    readname = 'bogus-genome-chr1_115_449_0:0:0_0:0:0_1f4/1'
    case = data_file('trio1/case1.fq')
    ctrls = kevlar.tests.data_glob('trio1/ctrl[1,2].fq')
    arglist = [
        'novel', '--ctrl-max', '0', '--case-min', '6', '--case', case,
        '--control', ctrls[0], '--control', ctrls[1], '--skip-until', readname
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.novel.main(args)
    out, err = capsys.readouterr()
    message = ('Found read bogus-genome-chr1_115_449_0:0:0_0:0:0_1f4/1 '
               '(skipped 1001 reads)')
    assert message in err
    assert '29 unique novel kmers in 14 reads' in err

    readname = 'BOGUSREADNAME'
    arglist = [
        'novel', '--ctrl-max', '0', '--case-min', '6', '--case', case,
        '--control', ctrls[0], '--control', ctrls[1], '--skip-until', readname
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.novel.main(args)
    out, err = capsys.readouterr()
    assert 'Found read' not in err
    assert '(skipped ' not in err
    assert 'Found 0 instances of 0 unique novel kmers in 0 reads' in err


def test_novel_output_has_mates():
    kid = data_file('microtrios/trio-na-proband.fq.gz')
    mom = data_file('microtrios/trio-na-mother.fq.gz')
    dad = data_file('microtrios/trio-na-father.fq.gz')
    testnovel = data_file('microtrios/novel-na.augfastq.gz')

    with NamedTemporaryFile(suffix='.augfastq') as novelfile:
        arglist = [
            'novel', '--out', novelfile.name, '--case', kid, '--case-min', '5',
            '--control', mom, '--control', dad, '--ctrl-max', '1',
            '--memory', '500K'
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.novel.main(args)

        intread_ids = set()
        mate_seqs = set()
        stream = kevlar.parse_augmented_fastx(kevlar.open(novelfile.name, 'r'))
        for read in stream:
            intread_ids.add(read.name)
            mate_seqs.update(read.mates)

        stream = kevlar.parse_augmented_fastx(kevlar.open(testnovel, 'r'))
        test_ids = set([r.name for r in stream])
        assert intread_ids == test_ids

        stream = kevlar.parse_augmented_fastx(kevlar.open(testnovel, 'r'))
        test_mate_seqs = set([m for r in stream for m in r.mates])
        assert mate_seqs == test_mate_seqs


def test_novel_save_counts():
    outdir = mkdtemp()
    try:
        for ind in ('father', 'mother', 'proband'):
            outfile = '{:s}/{:s}.ct'.format(outdir, ind)
            infile = data_file('microtrios/trio-na-{:s}.fq.gz'.format(ind))
            arglist = ['count', '--ksize', '27', '--memory', '500K', outfile,
                       infile]
            args = kevlar.cli.parser().parse_args(arglist)
            kevlar.count.main(args)

        arglist = [
            'novel', '--ksize', '27', '--out', outdir + '/novel.augfastq.gz',
            '--save-case-counts', outdir + '/kid.ct', '--save-ctrl-counts',
            outdir + '/mom.ct', outdir + '/dad.ct', '--case',
            data_file('microtrios/trio-na-proband.fq.gz'),
            '--control', data_file('microtrios/trio-na-mother.fq.gz'),
            '--control', data_file('microtrios/trio-na-father.fq.gz'),
            '--memory', '500K'
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.novel.main(args)

        counts = ('father', 'mother', 'proband')
        testcounts = ('dad', 'mom', 'kid')
        for c1, c2 in zip(counts, testcounts):
            f1 = '{:s}/{:s}.ct'.format(outdir, c1)
            f2 = '{:s}/{:s}.ct'.format(outdir, c2)
            assert filecmp.cmp(f1, f2)
    finally:
        rmtree(outdir)


def test_novel_save_counts_mismatch(capsys):
    outdir = mkdtemp()
    try:
        arglist = [
            'novel', '--ksize', '27', '--out', outdir + '/novel.augfastq.gz',
            '--save-case-counts', outdir + '/kid.ct', '--save-ctrl-counts',
            outdir + '/mom.ct', outdir + '/dad.ct', outdir + '/sibling.ct',
            '--case', data_file('microtrios/trio-k-proband.fq.gz'),
            '--control', data_file('microtrios/trio-k-mother.fq.gz'),
            '--control', data_file('microtrios/trio-k-father.fq.gz'),
            '--memory', '500K'
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.novel.main(args)
    finally:
        rmtree(outdir)

    out, err = capsys.readouterr()
    assert 'stubbornly refusing to save k-mer counts' in err


def test_novel_load_counts(capsys):
    file1 = data_file('simple-genome-case-reads.fa.gz')
    file2 = data_file('ambig.fasta')
    file3 = data_file('simple-genome-case.ct')
    file4, file5 = data_glob('simple-genome-ctrl?.ct')
    arglist = [
        'novel', '-k', '25',
        '--case', file1, file2, '--case-counts', file3,
        '--control-counts', file4, file5
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.novel.main(args)

    out, err = capsys.readouterr()
    assert 'counttables for 2 sample(s) provided' in err
