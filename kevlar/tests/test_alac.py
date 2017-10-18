#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import sys
import kevlar
from kevlar.tests import data_file


def test_pico_4(capsys):
    reads = data_file('pico-4.augfastq.gz')
    refr = data_file('human-random-pico.fa.gz')
    arglist = ['alac', '--ksize', '25', reads, refr]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.alac.main(args)
    out, err = capsys.readouterr()

    msg = 'assembled 22/28 reads from 1 connected component(s) into 1 contig'
    assert msg in err

    vcf = ('seq1\t1175768\t.\tT\tC\t.\tPASS\tVW=CCCTGCCATTATAGATGCTAGATTCACATC'
           'TTCATTTATTTTTACTTTT')
    assert vcf.strip() == out.strip()
