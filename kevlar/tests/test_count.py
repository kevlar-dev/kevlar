#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
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
from kevlar.tests import data_file, data_glob


@pytest.fixture
def triomask():
    mask = khmer.Counttable(19, 1e4, 4)
    mask.consume('TGAGGGGACTAGGTGATCAGGTGAGGGTTTCCCAGTTCCCGAAGATGACT')
    mask.consume('GATCTTTCGCTCCCTGTCATCAAGGAGTGATACGCGAAGTGCGTCCCCTT')
    mask.consume('GAAGTTTTGACAATTTACGTGAGCCCTACCTAGCGAAACAACAGAGAACC')
    return mask


@pytest.mark.parametrize('mask,numbands,band', [
    (None, None, None),
    (None, 9, 2),
    (triomask, None, None),
    (triomask, 23, 19),
])
def test_load_threading(mask, numbands, band):
    # Smoke test: make sure things don't explode when run in "threaded" mode.
    infiles = data_glob('trio1/case1.fq')
    sketch = kevlar.counting.load_sample_seqfile(
        infiles, 19, 1e7, mask=mask, numbands=numbands, band=band, numthreads=2
    )


def test_load_sketches():
    infiles = data_glob('test.counttable')
    sketches = kevlar.counting.load_samples_sketchfiles(infiles, maxfpr=0.5)
    for sketch in sketches:
        assert sketch.get('CCTGATATCCGGAATCTTAGC') > 0
        assert sketch.get('GATTACA' * 3) == 0


@pytest.mark.parametrize('numbands,band,kmers_stored', [
    (0, 0, 947),
    (2, 1, 500),
    (16, 7, 68),
])
def test_count_simple(numbands, band, kmers_stored, capsys):
    with NamedTemporaryFile(suffix='.counttable') as ctrl1out, \
            NamedTemporaryFile(suffix='.counttable') as ctrl2out, \
            NamedTemporaryFile(suffix='.counttable') as caseout:
        case = data_file('simple-genome-case-reads.fa.gz')
        ctrls = data_glob('simple-genome-ctrl[1,2]-reads.fa.gz')
        arglist = [
            'count',
            '--case', caseout.name, case,
            '--control', ctrl1out.name, ctrls[0],
            '--control', ctrl2out.name, ctrls[1],
            '--ksize', '25', '--memory', '5K', '--ctrl-max', '0',
            '--num-bands', str(numbands), '--band', str(band),
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.count.main(args)
    out, err = capsys.readouterr()

    assert '600 reads processed' in str(err)
    assert '{:d} distinct k-mers stored'.format(kmers_stored) in str(err)


def test_count_threading():
    with NamedTemporaryFile(suffix='.counttable') as ctrl1out, \
            NamedTemporaryFile(suffix='.counttable') as ctrl2out, \
            NamedTemporaryFile(suffix='.counttable') as caseout:
        case = data_file('trio1/case1.fq')
        ctrls = data_glob('trio1/ctrl[1,2].fq')
        arglist = [
            'count',
            '--ksize', '19', '--memory', '500K', '--threads', '2',
            '--case', caseout.name, case,
            '--control', ctrl1out.name, ctrls[0],
            '--control', ctrl2out.name, ctrls[1],
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.count.main(args)

    # No checks, just doing a "smoke test" to make sure things don't explode
    # when counting is done in "threaded" mode.
