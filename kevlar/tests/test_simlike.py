#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import khmer
import kevlar
from kevlar.tests import data_file
from kevlar.simlike import get_abundances


def test_get_abundances():
    kid = khmer.Counttable(31, 1e6, 4)
    mom = khmer.Counttable(31, 1e6, 4)
    dad = khmer.Counttable(31, 1e6, 4)
    ref = khmer.Counttable(31, 1e6, 4)
    kid.consume_seqfile(data_file('minitrio/trio-proband.fq.gz'))
    mom.consume_seqfile(data_file('minitrio/trio-mother.fq.gz'))
    dad.consume_seqfile(data_file('minitrio/trio-father.fq.gz'))
    ref.consume_seqfile(data_file('minitrio/refr.fa'))

    window = 'ACTACCCTAACTTTGGAAGACTATCAAAAACCCATTTCTGGGGGTGGAGGGGAGGGAGACA'
    abund, ndropped = get_abundances(window, kid, (mom, dad), ref)
    assert ndropped == 0
    assert abund == [
        [7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 7, 7, 9, 9, 9, 9, 8, 8, 9, 7, 6, 6, 6,
         6, 6, 6, 6, 6, 6, 6, 7],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
         1, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1, 0, 0, 0, 0, 0, 0, 0],
    ]
