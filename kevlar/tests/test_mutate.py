#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from sys import stderr
import pytest
import kevlar
from kevlar.mutate import Mutation
from kevlar.tests import data_file


def test_load_mutations_x():
    instream = kevlar.open(data_file('muts-x.txt'), 'r')
    mutations = kevlar.mutate.load_mutations(instream, stderr)
    assert len(mutations) == 1
    assert '1' in mutations
    assert len(mutations['1']) == 1

    mut = mutations['1'][0]
    assert mut == Mutation(seq='1', pos=441274, type='snv', data='3')
