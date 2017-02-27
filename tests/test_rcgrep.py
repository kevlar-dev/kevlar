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

    args.query = ['ACCTTATTAAGTCACGCCC']
    args.file = ['tests/data/collect.beta.1.txt']
    kevlar.rcgrep.main(args)

    out, err = capsys.readouterr()
    for line in out.split('\t'):
        print(line)
        assert args.query[0] in line or kevlar.revcom(args.query[0]) in line
