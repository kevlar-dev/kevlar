#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from io import StringIO
from shutil import rmtree, copyfile
import sys
from tempfile import mkdtemp

import kevlar
from kevlar.reference import autoindex, load_refr_cutouts, ReferenceCutout
from kevlar.reference import KevlarBWAError, KevlarInvalidCutoutDeflineError
from kevlar.reference import KevlarDeflineSequenceLengthMismatchError
from kevlar.tests import data_file
import pytest


def test_cutout_basic():
    c1 = ReferenceCutout()
    assert c1.interval == (None, None, None)

    c2 = ReferenceCutout('1_1000-2000')
    assert c2.defline == '1_1000-2000'
    assert c2.sequence is None
    assert c2.interval == ('1', 1000, 2000)

    with pytest.raises(KevlarInvalidCutoutDeflineError):
        c3 = ReferenceCutout('deFlIne FOrMaT WHat arEYoutALKingAb out')

    c4 = ReferenceCutout('chr3_1000-2000', 'A' * 1000)
    assert c4.defline == 'chr3_1000-2000'
    assert c4.sequence == 'A' * 1000

    assert c4.local_to_global(40) == 1040

    with pytest.raises(KevlarDeflineSequenceLengthMismatchError):
        c5 = ReferenceCutout('scaffold_4000-5000', 'A' * 42)


def test_load_cutouts():
    instream = kevlar.open(data_file('ssc218.gdna.fa'), 'r')
    cutouts = list(load_refr_cutouts(instream))
    assert len(cutouts) == 1
    assert cutouts[0].defline == '6_23229978-23230336'
    assert cutouts[0].sequence.startswith('GAACTCTCAATAAGGAATGTAATTAGAGTCATGT')
    assert cutouts[0].sequence.endswith('GTTAAACAATGGATACAAAATTGATAGAAACAATTA')


def test_autoindex():
    # Bogus directory
    with pytest.raises(KevlarBWAError) as e:
        autoindex('/a/truly/bogus/dir/seqs.fa')
    assert 'does not exist' in str(e)

    tmpdir = mkdtemp()
    try:
        # Real directory, non-existent file
        filename = tmpdir + '/refr.fa'
        with pytest.raises(KevlarBWAError) as e:
            autoindex(filename)
        assert 'does not exist' in str(e)

        # Should successfully index the sequence
        log = StringIO()
        copyfile(data_file('bogus-genome/refr.fa'), filename)
        autoindex(filename, logstream=log)
        assert 'BWA index not found' in log.getvalue()
        assert 'indexing now' in log.getvalue()

        # Should find existing index
        log = StringIO()
        autoindex(filename, logstream=log)
        assert log.getvalue() == ''
    finally:
        rmtree(tmpdir)


def test_bwa_align_coords():
    seq = ('TGACGTGACCCCAAGAAAACACTGCACCCAACTTCTTTCTTTAAGCCTTCGTGTGTGCAGGAGGAG'
           'GCCAGCCCTGGTTTCAAAATTGTTCCTCAGCATT')
    fasta = '>seq1\n{}\n'.format(seq)
    args = ['bwa', 'mem', data_file('fiveparts-refr.fa.gz'), '-']
    aligner = kevlar.reference.bwa_align(args, fasta)
    mappings = list(aligner)
    assert len(mappings) == 1
    assert mappings[0] == ('seq1', 50, 150, seq)


def test_bwa_failure():
    args = ['bwa', 'mem', data_file('not-a-real-file.fa'), '-']
    with pytest.raises(KevlarBWAError) as e:
        aligner = kevlar.reference.bwa_align(args, '>seq1\nACGT')
        pos = list(aligner)
