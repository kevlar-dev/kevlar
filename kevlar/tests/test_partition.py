#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import os
import pytest
import re
import shutil
import sys
import tempfile
import kevlar
from kevlar.partition import partition


def test_mkdirp():
    tempdir = tempfile.mkdtemp()

    path1 = os.path.join(tempdir, 'partitions')
    assert kevlar.mkdirp(path1) == path1

    path2 = os.path.join(tempdir, 'partitions2', 'part')
    path2test = os.path.join(tempdir, 'partitions2')
    assert kevlar.mkdirp(path2, trim=True) == path2test

    path3 = os.path.join(tempdir, 'partitions3', 'a', 'long', 'path', 'PART')
    path3test = os.path.join(tempdir, 'partitions3', 'a', 'long', 'path')
    assert kevlar.mkdirp(path3, trim=True) == path3test

    shutil.rmtree(tempdir)


def test_partition_dedup(capsys):
    infile = kevlar.tests.data_file('dup.augfastq')
    tempdir = tempfile.mkdtemp()

    arglist = ['partition', '--split', tempdir + '/dedup', infile]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.partition.main(args)
    out, err = capsys.readouterr()
    assert 'grouped 16 reads into 1 connected components' in err

    outfile = tempdir + '/dedup.cc1.augfastq.gz'
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

    arglist = ['partition', '--no-dedup', '--split', tempdir + '/nodedup',
               infile]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.partition.main(args)
    out, err = capsys.readouterr()
    assert 'grouped 18 reads into 1 connected components' in err

    outfile = tempdir + '/nodedup.cc1.augfastq.gz'
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


def test_partition_dedup_minabund(capsys):
    infile = kevlar.tests.data_file('dupl-part.augfastq.gz')
    readstream = kevlar.parse_augmented_fastx(kevlar.open(infile, 'r'))
    partitioner = partition(readstream, minabund=5)
    partitions = list(partitioner)
    assert len(partitions) == 0


def test_partition_nodedup_minabund(capsys):
    infile = kevlar.tests.data_file('dupl-part-2reads.augfastq.gz')
    readstream = kevlar.parse_augmented_fastx(kevlar.open(infile, 'r'))
    partitioner = partition(readstream, minabund=5, dedup=False)
    partitions = list(partitioner)
    assert len(partitions) == 0


def test_pico_ccs():
    infile = kevlar.tests.data_file('pico-filtered.fq.gz')
    readstream = kevlar.parse_augmented_fastx(kevlar.open(infile, 'r'))
    partitions = [p for p in partition(readstream, minabund=6)]
    assert len(partitions) == 10


def test_pico_stream(capsys):
    infile = kevlar.tests.data_file('pico-filtered.fq.gz')
    arglist = ['partition', '--min-abund', '5', infile]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.partition.main(args)
    out, err = capsys.readouterr()

    kvcc_labels = set()
    for line in out.strip().split('\n'):
        kvcc_match = re.search(r'(kvcc=\d+)', line)
        if kvcc_match:
            kvcc_labels.add(kvcc_match.group(1))
    assert len(kvcc_labels) == 10
