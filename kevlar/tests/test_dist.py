#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import filecmp
import sys
from tempfile import NamedTemporaryFile

import pytest
import kevlar
from kevlar.dist import count_first_pass
from kevlar.tests import data_file
from khmer import Counttable, Nodetable


def test_count_first_pass():
    mask = Nodetable.load(data_file('minitrio/mask.nt'))
    counts = Counttable(31, 1e4, 4)
    seqfile = data_file('minitrio/trio-proband.fq.gz')
    count_first_pass([seqfile], counts, mask)
    with NamedTemporaryFile(suffix='.ct') as countfile:
        counts.save(countfile.name)
        testcountfile = data_file('minitrio/trio-proband-mask-counts.ct')
        assert filecmp.cmp(testcountfile, countfile.name) is True
