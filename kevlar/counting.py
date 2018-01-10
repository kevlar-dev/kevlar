#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import khmer
import kevlar


def load_samples(samplelists, ksize, memory, mask=None, memfraction=None,
                 maxfpr=0.2, maxabund=1, numbands=None, band=None,
                 outfiles=None, numthreads=1, logfile=sys.stderr):
    """
    Load a group of related samples using a memory-efficient strategy.

    By default, each sample is loaded into a dedicated counttable, which occupy
    `memory` bytes of memory each. Setting `memfraction` to a value between 0.0
    and 1.0 will activate "masked" mode.

    If `mask` is provided, it serves as a mask for all other samples. If it is
    not provided, the first sample is loaded normally and then serves as a mask
    for all subsequent samples.

    In "masked mode", sample uses only `memory * memfraction` bytes of memory,
    and any k-mer present in the mask (above a given threshold `maxabund`) is
    ignored. In this way, we avoid taking up space storing abundances for
    k-mers we know we're not interested in.
    """
    numsamples = len(samplelists)
    if outfiles is None:
        outfiles = [None] * numsamples
    if numsamples != len(outfiles):
        message = '# of samples ({:d}) '.format(numsamples)
        message += 'does not match # of outfiles ({:d})'.format(len(outfiles))
        raise KevlarOutfileMismatchError(message)
    message = 'computing k-mer abundances for {:d} samples'.format(numsamples)
    print('[kevlar::counting]    ', message, file=logfile)

    sketches = list()
    for seqfiles, outfile in zip(samplelists, outfiles):
        sketchmem = memory if memfraction is None else memory * memfraction
        mymask = mask
        if memfraction is not None and len(sketches) == 0 and mask is None:
            mymask = sketches[0]
        sketch = kevlar.count.load_sample_seqfile(
            seqfiles, ksize, sketchmem, maxfpr=maxfpr, mask=mymask,
            numbands=numbands, band=band, outfile=outfile,
            numthreads=numthreads, logfile=logfile
        )
        sketches.append(sketch)
    return sketches


def load_samples_sketchfiles(sketchfiles, maxfpr=0.2, logfile=sys.stderr):
    """Load samples from pre-computed k-mer abundances."""
    sketches = list()
    for sketchfile in sketchfiles:
        message = 'loading sketchfile "{}"...'.format(sketchfile)
        print('[kevlar::counting]    ', message, end='', file=logfile)
        sketch = kevlar.sketch.autoload(sketchfile)
        fpr = kevlar.sketch.estimate_fpr(sketch)
        message = 'done! estimated false positive rate is {:1.3f}'.format(fpr)
        if fpr > maxfpr:
            message += ' (FPR too high, bailing out!!!)'
            raise SystemExit(message)
        print(message, file=logfile)
        sketches.append(sketch)
    return sketches
