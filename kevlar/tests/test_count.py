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
import sys


@pytest.mark.parametrize('numbands,band,kmers_stored', [
    (0, 0, 15600),
    (2, 1, 7663),
    (16, 7, 859),
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
            '--ksize', '25',  '--ctrl-max', '0',
            '--mem-frac', '0.1', '--memory', '50K',
            '--num-bands', str(numbands), '--band', str(band),
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.count.main(args)
    out, err = capsys.readouterr()

    print(err, file=sys.stderr)
    assert '600 reads processed' in str(err)
    assert '{:d} k-mers stored'.format(kmers_stored) in str(err)
