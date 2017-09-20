#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict

import khmer
from khmer import khmer_args
import kevlar


def split_infiles_outfiles(filelist):
    outfiles = [flist[0] for flist in filelist]
    infilelists = [flist[1:] for flist in filelist]
    return outfiles, infilelists


def main(args):
    if (args.num_bands is None) is not (args.band is None):
        raise ValueError('Must specify --num-bands and --band together')
    myband = args.band - 1 if args.band else None

    timer = kevlar.Timer()
    timer.start()

    timer.start('loadctrl')
    print('[kevlar::count] Loading control samples', file=args.logfile)
    outfiles, infilelists = split_infiles_outfiles(args.control)
    controls = kevlar.counting.load_samples(
        infilelists, args.ksize, args.memory, outfiles=outfiles,
        memfraction=args.mem_frac, maxfpr=args.max_fpr, maxabund=args.ctrl_max,
        mask=None, numbands=args.num_bands, band=myband,
        numthreads=args.threads, logfile=args.logfile
    )
    elapsed = timer.stop('loadctrl')
    numcontrols = len(controls)
    message = '{:d} samples loaded in {:.2f} sec'.format(numcontrols, elapsed)
    print('[kevlar::count]', message, file=args.logfile)

    print('[kevlar::count] Loading case samples', file=args.logfile)
    timer.start('loadcase')
    outfiles, infilelists = split_infiles_outfiles(args.case)
    casemask = outfiles[0] if args.mem_frac else None
    cases = kevlar.counting.load_samples(
        infilelists, args.ksize, args.memory, outfiles=outfiles,
        memfraction=args.mem_frac, maxfpr=args.max_fpr, maxabund=args.ctrl_max,
        mask=casemask, numbands=args.num_bands, band=myband,
        numthreads=args.threads, logfile=args.logfile
    )
    elapsed = timer.stop('loadcase')
    numcases = len(cases)
    message = '{:d} sample(s) loaded in {:.2f} sec'.format(numcases, elapsed)
    print('[kevlar::count]', message, file=args.logfile)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::count]', message, file=args.logfile)
