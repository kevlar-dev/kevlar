#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import glob
import pytest
import re
import screed
import kevlar
from kevlar.tests import data_file, data_glob


@pytest.mark.parametrize('numbands,band,kmers_stored', [
    (0, 0, 15600),
    (2, 1, 7992),
    (16, 7, 1218),
])
def test_count_simple(numbands, band, kmers_stored, capsys):
    case = data_file('simple-genome-case-reads.fa.gz')
    ctrls = data_glob('simple-genome-ctrl[1,2]-reads.fa.gz')
    arglist = ['count', '--ksize', '25', '--memory', '5K', '--ctrl_max', '0',
               '--num-bands', str(numbands), '--band', str(band),
               '--case', case, '--controls'] + ctrls
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.count.main(args)
    out, err = capsys.readouterr()

    assert '600 reads processed' in str(err)
    assert '{:d} k-mers stored'.format(kmers_stored) in str(err)
