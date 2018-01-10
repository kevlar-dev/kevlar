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
import kevlar


def main(args):
    if (args.num_bands is None) is not (args.band is None):
        raise ValueError('Must specify --num-bands and --band together')
    myband = args.band - 1 if args.band else None

    timer = kevlar.Timer()
    timer.start()

    sketch = kevlar.counting.load_sample_seqfile(
        args.seqfile, args.ksize, args.memory, args.max_fpr,
        numbands=args.num_bands, band=myband, numthreads=args.threads,
        outfile=args.counttable, logfile=args.logfile
    )

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::count]', message, file=args.logfile)
