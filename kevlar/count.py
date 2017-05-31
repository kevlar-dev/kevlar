#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict

import khmer
from khmer import khmer_args
import kevlar


def main(args):
    if (args.num_bands is None) is not (args.band is None):
        raise ValueError('Must specify --num-bands and --band together')

    timer = kevlar.Timer()
    timer.start()

    timer.start('loadctrl')
    print('[kevlar::count] Loading control samples', file=args.logfile)
    controls = kevlar.counting.load_samples_with_dilution(
        args.control, args.ksize, args.memory, memfraction=args.mem_frac,
        maxfpr=args.max_fpr, maxabund=args.ctrl_max, masks=None,
        numbands=args.num_bands, band=args.band, logfile=args.logfile
    )
    elapsed = timer.stop('loadctrl')
    numcontrols = len(controls)
    message = '{:d} samples loaded in {:.2f} sec'.format(numcontrols, elapsed)
    print('[kevlar::count]', message, file=args.logfile)

    print('[kevlar::count] Loading case samples', file=args.logfile)
    timer.start('loadcase')
    cases = kevlar.counting.load_samples_with_dilution(
        args.case, args.ksize, args.memory, memfraction=args.mem_frac,
        maxfpr=args.max_fpr, maxabund=args.ctrl_max, masks=controls,
        numbands=args.num_bands, band=args.band, logfile=args.logfile
    )
    elapsed = timer.stop('loadcase')
    numcases = len(cases)
    message = '{:d} sample(s) loaded in {:.2f} sec'.format(numcases, elapsed)
    print('[kevlar::count]', message, file=args.logfile)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::count]', message, file=args.logfile)
