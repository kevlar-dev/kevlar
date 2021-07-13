#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import kevlar
from kevlar import ProgressIndicator
from kevlar.tests import data_file
import pytest
from time import sleep


def test_progress(capsys):
    kevlar.logstream, logstream = sys.stderr, kevlar.logstream
    message = 'processed {counter} partitions'
    logger = ProgressIndicator(message, interval=1, breaks=[10, 100, 1000])
    for _ in range(12000):
        logger.update()
    kevlar.logstream = logstream
    out, err = capsys.readouterr()
    print(err)
    assert open(data_file('progind.txt'), 'r').read().strip() == err.strip()


def test_progress_timer(capsys):
    kevlar.logstream, logstream = sys.stderr, kevlar.logstream
    message = 'loaded {counter} reads'
    logger = ProgressIndicator(message, interval=1, breaks=[10], usetimer=True)
    for _ in range(111):
        sleep(0.01)
        logger.update()
    kevlar.logstream = logstream
    out, err = capsys.readouterr()
    print(err)
    loglines = err.strip().split('\n')
    timelines = [ln for ln in loglines if 'seconds elapsed' in ln]
    assert len(loglines) == 20
    assert len(timelines) == 20
