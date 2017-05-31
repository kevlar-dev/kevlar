#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import kevlar
from kevlar.tests import data_file


def test_basic(capsys):
    args = kevlar.cli.parser().parse_args([
        'reaugment',
        data_file('reaugment.augfastq'),
        data_file('reaugment.fq')
    ])
    kevlar.reaugment.main(args)

    out, err = capsys.readouterr()
    testout = open(data_file('reaugment.out'), 'r').read()
    assert out == testout
