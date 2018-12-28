#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 Battelle National Biodefense Institute
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
from kevlar.tests import data_glob
import pytest


def test_unband_beta():
    infiles = data_glob('collect.beta.?.txt')
    instream = kevlar.seqio.afxstream(infiles)
    merger = kevlar.unband.unband(instream, numbatches=2)
    reads = list(merger)
    reads.sort(key=lambda r: r.name)
    assert len(reads) == 8
    assert len(reads[0].annotations) == 4


def test_unband_helium():
    infiles = data_glob('helium-unband/novel.band?.augfastq.gz')
    instream = kevlar.seqio.afxstream(infiles)
    merger = kevlar.unband.unband(instream, numbatches=16)
    reads = list(merger)
    assert len(reads) == 135

    readname = 'seq1_haplo1_285110_285519_1:0:0_0:0:0_2dbcd/1'
    someread = [r for r in reads if r.name == readname][0]
    assert len(someread.annotations) == 75


def test_unband_cli(capsys):
    infiles = data_glob('helium-unband/novel.band?.augfastq.gz')
    arglist = ['unband'] + infiles
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.unband.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    qualdeflines = [ln for ln in outlines if ln == '+']
    assert len(qualdeflines) == 135
