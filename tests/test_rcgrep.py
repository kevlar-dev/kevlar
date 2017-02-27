#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import sys

import khmer
import kevlar


def test_simple(capsys):
    args = type('', (), {})()

    args.before = 0
    args.after = 0
    args.query = ['ACCTTATTAAGTCACGCCC']
    args.file = ['tests/data/collect.beta.1.txt']
    kevlar.rcgrep.main(args)

    out, err = capsys.readouterr()
    for line in out.split('\t'):
        assert args.query[0] in line or kevlar.revcom(args.query[0]) in line


def test_before(capsys):
    args = type('', (), {})()

    args.before = 1
    args.after = 0
    args.query = ['GTATACTACTGCGGCATGGGA']
    args.file = ['tests/data/trio1/case2.fq']
    kevlar.rcgrep.main(args)

    out, err = capsys.readouterr()
    outlines = out.split('\t')
    for line1, line2 in zip(outlines[::2], outlines[1::2]):
        assert 'chr1' in line1
        assert args.query[0] in line2 or kevlar.revcom(args.query[0]) in line2


def test_after(capsys):
    args = type('', (), {})()

    args.before = 0
    args.after = 1
    args.query = ['GTATACTACTGCGGCATGGGA']
    args.file = ['tests/data/trio1/case2.fq']
    kevlar.rcgrep.main(args)

    out, err = capsys.readouterr()
    outlines = out.split('\t')
    for line1, line2 in zip(outlines[::2], outlines[1::2]):
        assert args.query[0] in line1 or kevlar.revcom(args.query[0]) in line1
        assert '+' == line2


def test_before_and_after(capsys):
    args = type('', (), {})()

    args.before = 1
    args.after = 1
    args.query = ['GTATACTACTGCGGCATGGGA']
    args.file = ['tests/data/trio1/case2.fq']
    kevlar.rcgrep.main(args)

    out, err = capsys.readouterr()
    outlines = out.split('\t')
    for l1, l2, l3 in zip(outlines[::2], outlines[1::2], outlines[2::2]):
        assert 'chr1' in l1
        assert args.query[0] in l2 or kevlar.revcom(args.query[0]) in l2
        assert '+' == l3
