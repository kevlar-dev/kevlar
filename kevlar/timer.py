#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import time


class Timer(object):
    def __init__(self):
        self._start_times = dict()
        self._stop_times = dict()

    def start(self, key=None):
        if key is None:
            key = ''
        if key in self._start_times:
            raise ValueError('Timer already started for "' + key + '"')
        self._start_times[key] = time.time()

    def stop(self, key=None):
        if key is None:
            key = ''
        if key not in self._start_times:
            raise ValueError('No timer started for "' + key + '"')
        self._stop_times[key] = time.time()
        return self._stop_times[key] - self._start_times[key]

    def probe(self, key=None):
        if key is None:
            key = ''
        if key not in self._start_times:
            raise ValueError('No timer started for "' + key + '"')
        current = time.time()
        return current - self._start_times[key]
