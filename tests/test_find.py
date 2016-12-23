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
import sys
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import kevlar


@pytest.fixture
def trio_args():
    args = type('', (), {})()

    args.controls = glob.glob('tests/data/trio1/ctrl[1,2].fq')
    args.ctrl_max = 0
    args.case_min = 8
    args.ksize = 13
    args.graph_memory = 1e6
    args.out = StringIO()
    args.flush = False
    args.kmers_out = None
    args.paths_out = StringIO()
    args.collapse = True
    args.upint = 1000
    args.logfile = StringIO()
    args.case = 'tests/data/trio1/case1.fq'

    return args


def run_single(case, ctrl, gmem, mutseq, ksize, trio_args):
    trio_args.ksize = ksize
    trio_args.graph_memory = gmem
    trio_args.case = 'tests/data/trio1/{}.fq'.format(case)
    trio_args.controls = glob.glob('tests/data/trio1/{}.fq'.format(ctrl))
    kevlar.find.main(trio_args)

    pathdata = trio_args.paths_out.getvalue()
    paths = sorted([ln.split(',')[0] for ln in pathdata.strip().split('\n')])
    assert len(paths) == 1

    mutseqrc = kevlar.revcom(mutseq)
    assert mutseq in paths[0] or mutseqrc in paths[0]


@pytest.mark.parametrize('case,ctrl,gmem,mutseq,ksize', [
    ('case1', 'ctrl[1,2]', 5e5, 'CCGCACCATTT', 13),
    ('case1', 'ctrl[1,2]', 1e6, 'CCGCACCATTT', 13),
    ('case4', 'ctrl[1,2]', 5e5, 'GGTCAATAGG', 13),
    ('case4', 'ctrl[1,2]', 1e6, 'GGTCAATAGG', 13),
    ('case5', 'ctrl[3,4]', 1e6, 'GGTCAATAGG', 13),
    ('case6', 'ctrl[5,6]', 1e6, 'GTCAATAG', 13),
])
def test_find_single_mutation(case, ctrl, gmem, mutseq, ksize, trio_args):
    run_single(case, ctrl, gmem, mutseq, ksize, trio_args)


@pytest.mark.long
@pytest.mark.parametrize('case,ctrl,gmem,mutseq,ksize', [
    ('case1', 'ctrl[1,2]', 1e6, 'CCGCACCATTT', 11),
    ('case1', 'ctrl[1,2]', 1e6, 'CCGCACCATTT', 15),
    ('case1', 'ctrl[1,2]', 1e6, 'CCGCACCATTT', 17),
    ('case1', 'ctrl[1,2]', 1e6, 'CCGCACCATTT', 19),
    ('case1', 'ctrl[1,2]', 5e5, 'CCGCACCATTT', 11),
    ('case1', 'ctrl[1,2]', 5e5, 'CCGCACCATTT', 15),
    ('case1', 'ctrl[1,2]', 5e5, 'CCGCACCATTT', 17),
    ('case1', 'ctrl[1,2]', 5e5, 'CCGCACCATTT', 19),
    ('case4', 'ctrl[1,2]', 5e5, 'GGTCAATAGG', 11),
    ('case4', 'ctrl[1,2]', 5e5, 'GGTCAATAGG', 15),
    ('case4', 'ctrl[1,2]', 5e5, 'GGTCAATAGG', 17),
    ('case4', 'ctrl[1,2]', 5e5, 'GGTCAATAGG', 19),
    ('case4', 'ctrl[1,2]', 1e6, 'GGTCAATAGG', 11),
    ('case4', 'ctrl[1,2]', 1e6, 'GGTCAATAGG', 15),
    ('case4', 'ctrl[1,2]', 1e6, 'GGTCAATAGG', 17),
    ('case4', 'ctrl[1,2]', 1e6, 'GGTCAATAGG', 19),
    ('case5', 'ctrl[3,4]', 1e6, 'GGTCAATAGG', 11),
    ('case5', 'ctrl[3,4]', 1e6, 'GGTCAATAGG', 15),
    ('case5', 'ctrl[3,4]', 1e6, 'GGTCAATAGG', 17),
    ('case5', 'ctrl[3,4]', 1e6, 'GGTCAATAGG', 19),
])
def test_find_single_mutation_long(case, ctrl, gmem, mutseq, ksize, trio_args):
    run_single(case, ctrl, gmem, mutseq, ksize, trio_args)


def test_find_case2(trio_args):
    trio_args.case = 'tests/data/trio1/case2.fq'
    trio_args.ksize = 13
    kevlar.find.main(trio_args)

    pathdata = trio_args.paths_out.getvalue()
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


def test_find_case3(trio_args):
    trio_args.case = 'tests/data/trio1/case3.fq'
    kevlar.find.main(trio_args)

    pathdata = trio_args.paths_out.getvalue()
    pathlines = [ln for ln in pathdata.strip().split('\n') if 'chr1' in ln]
    paths = sorted([ln.split(',')[0] for ln in pathlines])
    assert len(paths) == 1

    #              ┌----------- SNV here
    mutseq = 'CCGCACCATTT'
    mutseqrc = kevlar.revcom(mutseq)
    assert mutseq in paths[0] or mutseqrc in paths[0]
