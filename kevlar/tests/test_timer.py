#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import time
import pytest
import kevlar
from kevlar.timer import Timer


def test_testy_mctestface(capsys):
    t = Timer()
    t.start()
    t.start('task1')
    time.sleep(0.1)
    assert t.probe() > 0
    assert t.probe('task1') > 0.0
    with pytest.raises(ValueError) as ve:
        t.probe('task2')
    assert 'No timer started for "task2"' in str(ve)
    elapsed = t.stop('task1')
    assert elapsed > 0.0
    t.start('task2')
    time.sleep(0.1)
    assert t.probe('task2') > 0.0
    elapsed = t.stop('task2')
    assert elapsed > 0.0
    with pytest.raises(ValueError) as ve:
        t.stop('task3')
    assert 'No timer started for "task3"' in str(ve)
    t.start('task3')
    with pytest.raises(ValueError) as ve:
        t.start('task3')
    assert 'Timer already started for "task3"' in str(ve)
