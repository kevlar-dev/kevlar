#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import filecmp
from os import remove
import pytest
import subprocess
import sys
from tempfile import NamedTemporaryFile, mkdtemp
from shutil import rmtree
import kevlar
from kevlar.tests import data_file, data_glob


@pytest.fixture
def pico_trio(request):
    mother = NamedTemporaryFile(suffix='.fq.gz', delete=False)
    father = NamedTemporaryFile(suffix='.fq.gz', delete=False)
    proband = NamedTemporaryFile(suffix='.fq.gz', delete=False)
    refr = NamedTemporaryFile(suffix='.fq.gz', delete=False)
    filenames = [mother.name, father.name, proband.name, refr.name]

    def teardown():
        for filename in filenames:
            remove(filename)
    request.addfinalizer(teardown)

    urls = [
        'https://osf.io/e8jb3/download?version=1',
        'https://osf.io/fuaty/download?version=1',
        'https://osf.io/f5trh/download?version=1',
        'https://osf.io/58rwa/download?version=1',
    ]
    for url, filename in zip(urls, filenames):
        subprocess.check_call(['curl', '-L', '-o', filename, url])
    subprocess.check_call(['ls', '-lhp'] + filenames)

    subprocess.check_call(['bwa', 'index', refr.name])
    return mother.name, father.name, proband.name, refr.name


@pytest.mark.toolong
def test_simplex_pico(pico_trio, capsys):
    mother, father, proband, refr = pico_trio

    with capsys.disabled():
        cases = kevlar.novel.load_samples(
            None, [[proband]], ksize=25, memory=5e6, maxfpr=0.6, numthreads=2,
            logstream=sys.stderr
        )
        controls = kevlar.novel.load_samples(
            None, [[mother], [father]], ksize=25, memory=5e6, maxfpr=0.6,
            numthreads=2, logstream=sys.stderr
        )
        mask = kevlar.filter.load_mask([refr], 25, 5e7, maxfpr=0.005,
                                       logstream=sys.stderr)

        caserecords = kevlar.multi_file_iter_screed([proband])
        workflow = kevlar.simplex.simplex(
            caserecords, cases[0], controls, refr, ksize=25, ctrlmax=0,
            casemin=6, mask=mask, filtermem=1e7, filterfpr=0.005,
            refrsctmem=1e8, logstream=sys.stderr
        )
        variants = [v for v in workflow]
    variants = sorted(variants, key=lambda v: v._pos)
    startpos = [v._pos + 1 for v in variants]
    teststartpos = [4073, 185752, 226611, 636699, 834646, 901124, 1175768,
                    1527139, 1631013, 2265795]
    assert len(variants) == 10
    assert startpos == teststartpos


def test_simplex_trio1(capsys):
    case = data_file('trio1/case1.fq')
    controls = data_glob('trio1/ctrl[1,2].fq')
    refr = data_file('bogus-genome/refr.fa')
    arglist = [
        'simplex', '--case', case, '--control', controls[0], '--control',
        controls[1], '--case-min', '6', '--ctrl-max', '0', '--novel-memory',
        '1M', '--novel-fpr', '0.2', '--filter-memory', '50K', '--mask-files',
        refr, '--mask-memory', '1M', '--filter-fpr', '0.005', '--ksize', '21',
        '--seed-size', '31', '--delta', '25',
        '--labels', 'Case', 'Ctrl1', 'Ctrl2', '--refr-sct-mem', '500M', refr
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.simplex.main(args)

    out, err = capsys.readouterr()
    # grep -v ^'#' out
    out = '\n'.join([l for l in out.split('\n') if not l.startswith('#')])

    testvcf = '\t'.join([
        'bogus-genome-chr1', '3567', '.', 'A', 'C', '.', 'PASS',
        'ALTWINDOW=GAAGGGCACACCTAACCGCACCATTTGCCGTGGAAGCATAA;CALLCLASS=1;'
        'CIGAR=25D95M25D;DROPPED=0;IKMERS=21;KSW2=82;LIKESCORE=252.625;'
        'LLDN=-62.504;LLFP=-778.456;LLIH=-315.129;PART=1;REFRWINDOW=GAAGGGCACA'
        'CCTAACCGCAACATTTGCCGTGGAAGCATAA;CONTIG=TTGGTGCCACGATCCGGCTATGGCGGAAGG'
        'GCACACCTAACCGCACCATTTGCCGTGGAAGCATAAAGGTCATCATTGAGGTGGTTCGTTCCGAT',
        'ALTABUND',
        '9,9,10,10,10,8,10,10,10,10,10,10,10,10,11,12,13,13,13,12,12',
        '0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0',
        '0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0'
    ])
    print(out.strip())
    print(testvcf)
    assert out.strip() == testvcf


def test_simplex_minitrio(capsys):
    proband = data_file('minitrio/trio-proband.fq.gz')
    mother = data_file('minitrio/trio-mother.fq.gz')
    father = data_file('minitrio/trio-father.fq.gz')
    refr = data_file('minitrio/refr.fa')

    arglist = [
        'simplex', '--novel-memory', '5M', '--case', proband,
        '--control', mother, '--control', father, '--threads', '2',
        '--ctrl-max', '1', '--case-min', '5', '--ksize', '25',
        refr
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.simplex.main(args)

    out, err = capsys.readouterr()
    print('DEBUG', out)
    # grep -v ^'#' out
    out = '\n'.join([l for l in out.split('\n') if not l.startswith('#')])
    assert len(out.strip().split('\n')) == 1


def test_simplex_save_counts():
    outdir = mkdtemp()
    try:
        for ind in ('father', 'mother', 'proband'):
            outfile = '{:s}/{:s}.ct'.format(outdir, ind)
            infile = data_file('microtrios/trio-k-{:s}.fq.gz'.format(ind))
            arglist = ['count', '--ksize', '27', '--memory', '500K', outfile,
                       infile]
            args = kevlar.cli.parser().parse_args(arglist)
            kevlar.count.main(args)

        arglist = [
            'simplex', '--ksize', '27', '--out', outdir + '/calls.vcf',
            '--save-case-counts', outdir + '/kid', '--save-ctrl-counts',
            outdir + '/mom', outdir + '/dad', '--case',
            data_file('microtrios/trio-k-proband.fq.gz'),
            '--control', data_file('microtrios/trio-k-mother.fq.gz'),
            '--control', data_file('microtrios/trio-k-father.fq.gz'),
            '--novel-memory', '500K', data_file('microtrios/refr-k.fa.gz')
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.simplex.main(args)

        counts = ('father', 'mother', 'proband')
        testcounts = ('dad', 'mom', 'kid')
        for c1, c2 in zip(counts, testcounts):
            f1 = '{:s}/{:s}.ct'.format(outdir, c1)
            f2 = '{:s}/{:s}.counttable'.format(outdir, c2)
            assert filecmp.cmp(f1, f2)
    finally:
        rmtree(outdir)
