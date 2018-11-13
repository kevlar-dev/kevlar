#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
import sys


class ProgressIndicator(object):
    """Progress indicator that prints log messages with decreasing frequency.

    Print progress updates frequently at first, but less frequently as the
    process continues. This gives a good idea of how quickly the processing is
    at first but ensures the progress updates don't overwhelm the terminal.
    """
    def __init__(self, message, interval=10, breaks=[100, 1000, 10000],
                 usetimer=False, logstream=sys.stderr):
        self.message = message
        self.counter = 0
        self.interval = interval
        self.nextupdate = interval
        self.breaks = breaks
        self.outstream = logstream
        self.timer = None
        if usetimer:
            self.timer = kevlar.Timer()
            self.timer.start()

    def update(self):
        if self.counter in self.breaks:
            self.interval = self.counter
        if self.counter >= self.nextupdate:
            self.nextupdate += self.interval
            message = self.message.format(counter=self.counter)
            if self.timer:
                elapsed = self.timer.probe()
                message += ' ({:.2f} seconds elapsed)'.format(elapsed)
            print(message, file=self.outstream)
        self.counter += 1
