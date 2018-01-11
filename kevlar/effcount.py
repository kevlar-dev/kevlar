#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import khmer
import kevlar


def load_samples(samplelists, ksize, memory, memfraction=None,
                 maxfpr=0.2, maxabund=1, numbands=None, band=None,
                 numthreads=1, logfile=sys.stderr):
    """
    Load a set of samples using a memory-efficient strategy.

    The first sample is loaded normally and occupies `memory` bytes of memory.
    This sample is then used as a mask for all subsequent samples, each of
    which occupies only `memory * memfraction` bytes of memory. Any k-mer
    absent from the mask (below the specified threshold `minabund`) is ignored
    in all subsequent samples. This allows one to avoid taking up space storing
    abundances for k-mers that are uninteresting.
    """
    numsamples = len(samplelists)
    assert numsamples > 1
    message = 'Computing k-mer abundances for {:d} samples'.format(numsamples)
    print('[kevlar::effcount]', message, file=logfile)

    bigsketch = kevlar.count.load_sample_seqfile(
        samplelists[0], ksize, memory, maxfpr=maxfpr, numbands=numbands,
        band=band, numthreads=numthreads, logfile=logfile
    )
    yield bigsketch

    for samplefilelist in samplelists[1:]:
        yield kevlar.count.load_sample_seqfile(
            samplefilelist, ksize, memory * memfraction, mask=bigsketch,
            maskmaxabund=maxabund, maxfpr=maxfpr, numbands=numbands, band=band,
            numthreads=numthreads, logfile=logfile
        )


def main(args):
    if (args.num_bands is None) is not (args.band is None):
        raise ValueError('Must specify --num-bands and --band together')
    myband = args.band - 1 if args.band else None
    if len(args.outfiles) != len(args.sample):
        message = 'number of outfiles must match number of declared samples'
        raise ValueError(message)

    timer = kevlar.Timer()
    timer.start()

    loader = load_samples(
        args.sample, args.ksize, args.memory, memfraction=args.memfrac,
        maxfpr=args.max_fpr, maxabund=args.max_abund, numbands=args.num_bands,
        band=args.band, numthreads=args.threads, logfile=args.logfile
    )
    for sketch, outfile in zip(loader, args.outfiles):
        sketch.save(outfile)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::effcount]', message, file=args.logfile)
