#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import shutil
import sys
import tempfile
import kevlar


def test_partition_dedup(capsys):
    infile = kevlar.tests.data_file('dup.augfastq')
    tempdir = tempfile.mkdtemp()

    arglist = ['partition', tempdir + '/dedup', infile]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.partition.main(args)
    out, err = capsys.readouterr()
    assert 'grouped 16 reads into 1 connected components' in err

    outfile = tempdir + '/dedup.cc0.augfastq.gz'
    stream = kevlar.open(outfile, 'r')
    parser = kevlar.parse_augmented_fastx(stream)
    readseqs = [r.sequence for r in parser]
    uniquereadseqs = set([kevlar.revcommin(s) for s in readseqs])
    testreads = [
        'AACGAACCACCTCAATGATGACCTTTATGCTTCCACGGCAAATGGTGCGG',
        'ACGAACCACCTCAATGATGACCTTTATGCTTCCACGGCAAATGGTGCGGT',
        'AGGGCACACCTAACCGCACCATTTGCCGTGGAAGCATAAAGGTCATCATT',
        'ATCGGAACGAACCACCTCAATGATGACCTTTATGCTTCCACGGCAAATGG',
        'CCACCTCAATGATGACCTTTATGCTTCCACGGCAAATGGTGCGGTTAGGT',
        'CCTCAATGATGACCTTTATGCTTCCACGGCAAATGGTGCGGTTAGGTGTG',
        'CGCACCATTTGCCGTGGAAGCATAAAGGTCATCATTGAGGTGGTTCGTTC',
        'CGGAAGGGCACACCTAACCGCACCATTTGCCGTGGAAGCATAAAGGTCAT',
        'CGGCTATGGCGGAAGGGCACACCTAACCGCACCATTTGCCGTGGAAGCAT',
        'CTATGGCGGAAGGGCACACCTAACCGCACCATTTGCCGTGGAAGCATAAA',
        'GCTTCCACGGCAAATGGTGCGGTTAGGTGTGCCCTTCCGCCATAGCCGGA',
        'GGAACGAACCACCTCAATGATGACCTTTATGCTTCCACGGCAAATGGTGC',
        'GGCAAATGGTGCGGTTAGGTGTGCCCTTCCGCCATAGCCGGATCGTGGCA',
        'TATGCTTCCACGGCAAATGGTGCGGTTAGGTGTGCCCTTCCGCCATAGCC',
        'TTATGCTTCCACGGCAAATGGTGCGGTTAGGTGTGCCCTTCCGCCATAGC',
        'TTGGTGCCACGATCCGGCTATGGCGGAAGGGCACACCTAACCGCACCATT',
    ]
    testreadseqs = set([kevlar.revcommin(s) for s in testreads])
    assert uniquereadseqs == testreadseqs

    shutil.rmtree(tempdir)


def test_partition_nodedup(capsys):
    infile = kevlar.tests.data_file('dup.augfastq')
    tempdir = tempfile.mkdtemp()

    arglist = ['partition', '--no-dedup', tempdir + '/nodedup', infile]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.partition.main(args)
    out, err = capsys.readouterr()
    assert 'grouped 18 reads into 1 connected components' in err

    outfile = tempdir + '/nodedup.cc0.augfastq.gz'
    stream = kevlar.open(outfile, 'r')
    parser = kevlar.parse_augmented_fastx(stream)
    readseqs = sorted([r.sequence for r in parser])
    assert readseqs == [
        'AACGAACCACCTCAATGATGACCTTTATGCTTCCACGGCAAATGGTGCGG',
        'ACGAACCACCTCAATGATGACCTTTATGCTTCCACGGCAAATGGTGCGGT',
        'AGGGCACACCTAACCGCACCATTTGCCGTGGAAGCATAAAGGTCATCATT',
        'ATCGGAACGAACCACCTCAATGATGACCTTTATGCTTCCACGGCAAATGG',
        'CACACCTAACCGCACCATTTGCCGTGGAAGCATAAAGGTCATCATTGAGG',
        'CCACCTCAATGATGACCTTTATGCTTCCACGGCAAATGGTGCGGTTAGGT',
        'CCTCAATGATGACCTTTATGCTTCCACGGCAAATGGTGCGGTTAGGTGTG',
        'CGCACCATTTGCCGTGGAAGCATAAAGGTCATCATTGAGGTGGTTCGTTC',
        'CGGAAGGGCACACCTAACCGCACCATTTGCCGTGGAAGCATAAAGGTCAT',
        'CGGCTATGGCGGAAGGGCACACCTAACCGCACCATTTGCCGTGGAAGCAT',
        'CTATGGCGGAAGGGCACACCTAACCGCACCATTTGCCGTGGAAGCATAAA',
        'GCTTCCACGGCAAATGGTGCGGTTAGGTGTGCCCTTCCGCCATAGCCGGA',
        'GGAACGAACCACCTCAATGATGACCTTTATGCTTCCACGGCAAATGGTGC',
        'GGCAAATGGTGCGGTTAGGTGTGCCCTTCCGCCATAGCCGGATCGTGGCA',
        'TATGCTTCCACGGCAAATGGTGCGGTTAGGTGTGCCCTTCCGCCATAGCC',
        'TTATGCTTCCACGGCAAATGGTGCGGTTAGGTGTGCCCTTCCGCCATAGC',
        'TTGGTGCCACGATCCGGCTATGGCGGAAGGGCACACCTAACCGCACCATT',
        'TTGGTGCCACGATCCGGCTATGGCGGAAGGGCACACCTAACCGCACCATT',
    ]

    shutil.rmtree(tempdir)
