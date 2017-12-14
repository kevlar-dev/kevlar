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
from kevlar import MutableString


@pytest.mark.parametrize('thestring', [
    ('GATTACA'),
    ('ATGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGA'),
    ('MAHERSHALALHASHBAZ')
])
def test_basic(thestring):
    ms = MutableString(thestring)
    assert str(ms) == thestring
    assert repr(ms) == thestring
    assert ms == thestring
    assert len(ms) == len(thestring)
    assert ms[3] == thestring[3]
    assert ms[2:6] == thestring[2:6]


def test_add_ops():
    ms = MutableString('HOWGOESIT')
    assert ms + 'MAN' == 'HOWGOESITMAN'
    assert ms + MutableString('MAN') == 'HOWGOESITMAN'

    ms += 'MANS!'
    assert ms == 'HOWGOESITMANS!'
    assert 'ITMANS' in ms

    ms[9:9] = "FELLOWHU"
    assert ms == 'HOWGOESITFELLOWHUMANS!'


def test_del_ops():
    ms = MutableString('Sesame seed buns')
    del ms[15]
    assert ms == 'Sesame seed bun'

    del ms[7:12]
    assert ms == 'Sesame bun'
