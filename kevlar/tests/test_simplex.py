#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from os import remove
import pytest
import subprocess
import sys
from tempfile import NamedTemporaryFile
import kevlar
from kevlar.tests import data_file


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


@pytest.mark.long
def test_simplex_pico(pico_trio):
    mother, father, proband, refr = pico_trio

    cases = kevlar.novel.load_samples(
        None, [[proband]], ksize=25, memory=5e6, maxfpr=0.6, numthreads=2,
    )
    controls = kevlar.novel.load_samples(
        None, [[mother], [father]], ksize=25, memory=5e6, maxfpr=0.6,
        numthreads=2,
    )
    mask = kevlar.filter.load_mask([refr], 25, 5e7, maxfpr=0.005)

    caserecords = kevlar.multi_file_iter_screed([proband])
    workflow = kevlar.simplex.simplex(
        caserecords, cases[0], controls, refr, ksize=25, ctrlmax=0, casemin=6,
        mask=mask, filtermem=1e7, filterfpr=0.005
    )
    variants = [v for v in workflow]
    variants = sorted(variants, key=lambda v: v._pos)
    startpos = [v._pos + 1 for v in variants]
    teststartpos = [4073, 185752, 226611, 636699, 834646, 901124, 1175768,
                    1527139, 1631013, 2265795]
    assert len(variants) == 10
    assert startpos == teststartpos
