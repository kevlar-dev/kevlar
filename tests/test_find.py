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
    args.controls = glob.glob('tests/data/trio1/ctrl[1,2].counts')
    args.ctrl_max = 0
    args.case_min = 8
    args.out = readout
    args.flush = False
    args.kmers_out = kmerout
    args.paths_out = pathout
    args.collapse = True
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

    pathdata = pathout.getvalue()
    paths = sorted([ln.split(',')[0] for ln in pathdata.strip().split('\n')])
    assert len(paths) == 1

    #              ┌----------- SNV here
    mutseq = 'CCGCACCATTT'
    mutseqrc = kevlar.revcom(mutseq)
    assert mutseq in paths[0] or mutseqrc in paths[0]


def test_find_case2(trio):
    args, pathout = trio
    args.case = 'tests/data/trio1/case2.counts'
    args.case_fastq = 'tests/data/trio1/case2.fq'
    kevlar.find.main(args)

    pathdata = pathout.getvalue()
    # This sorts by the value of the second field (which should sort by seq ID)
    pathlines = sorted(pathdata.strip().split('\n'),
                       key=lambda pd: pd.split('\t')[1])
    paths = [ln.split(',')[0] for ln in pathlines]
    assert len(paths) == 3

    mutseqs = [
        #      ┌----------- SNVs here
        'CCGCACCATTT',
        'AGAACTACCTT',
        'CATCCGTACGT',
    ]
    for mutseq, path in zip(mutseqs, paths):
        mutseqrc = kevlar.revcom(mutseq)
        assert mutseq in path or mutseqrc in path


def test_find_case3(trio):
    args, pathout = trio
    args.case = 'tests/data/trio1/case3.counts'
    args.case_fastq = 'tests/data/trio1/case3.fq'
    kevlar.find.main(args)

    pathdata = pathout.getvalue()
    pathlines = [ln for ln in pathdata.strip().split('\n') if 'chr1' in ln]
    paths = sorted([ln.split(',')[0] for ln in pathlines])
    assert len(paths) == 1

    #              ┌----------- SNV here
    mutseq = 'CCGCACCATTT'
    mutseqrc = kevlar.revcom(mutseq)
    assert mutseq in paths[0] or mutseqrc in paths[0]


def test_find_case4(trio):
    args, pathout = trio
    args.case = 'tests/data/trio1/case4.counts'
    args.case_fastq = 'tests/data/trio1/case4.fq'
    kevlar.find.main(args)

    pathdata = pathout.getvalue()
    paths = [ln.split(',')[0] for ln in pathdata.strip().split('\n')]
    assert len(paths) == 1

    #         ----┐┌--------- 5bp deletion between these nucleotides
    mutseq = 'GGTCAATAGG'
    mutseqrc = kevlar.revcom(mutseq)
    assert mutseq in paths[0] or mutseqrc in paths[0]
