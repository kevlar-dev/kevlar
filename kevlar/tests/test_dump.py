#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import sys
import kevlar
from kevlar.tests import data_file


def test_basic(capsys):
    arglist = [
        'dump',
        '--refr', data_file('bogus-genome/refr.fa'),
        data_file('bogus-genome/reads.bam'),
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.dump.main(args)
    out, err = capsys.readouterr()

    outputlines = out.strip().split('\n')
    assert len(outputlines) == 5 * 4  # 5 records, 4 lines per record
    assert 'read2' in outputlines[0]
    assert 'read4' in outputlines[4]
    assert 'read6' in outputlines[8]
    assert 'read7' in outputlines[12]
    assert 'read8' in outputlines[16]


def test_no_refr():
    bamstream = kevlar.open(data_file('bogus-genome/reads.bam'), 'r')
    records = [r for r in kevlar.dump.dump(bamstream)]
    assert len(records) == 8


def test_indels(capsys):
    arglist = [
        'dump',
        '-r', data_file('bogus-genome/refr.fa'),
        data_file('bogus-genome/reads-indels.bam'),
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.dump.main(args)
    out, err = capsys.readouterr()

    outputlines = out.strip().split('\n')
    assert len(outputlines) == 4 * 4  # 4 records, 4 lines per record


@pytest.mark.parametrize('strict,nrecords', [
    (True, 1),
    (False, 2)
])
def test_suffix(strict, nrecords):
    bamstream = kevlar.open(data_file('nopair.sam'), 'r')
    refrstream = kevlar.open(data_file('bogus-genome/refr.fa'), 'r')
    refr = kevlar.seqio.parse_seq_dict(refrstream)

    records = [r for r in kevlar.dump.dump(bamstream, refr, strict=strict)]
    assert len(records) == nrecords
    for record in records:
        assert record.name.endswith('/1') or record.name.endswith('/2')
