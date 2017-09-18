#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import re
import sys
import threading

import khmer
import kevlar


class KevlarSampleIOError(ValueError):
    pass


class KevlarOutfileMismatchError(ValueError):
    pass


def load_sample_seqfile(seqfiles, ksize, memory, maxfpr=0.2,
                        mask=None, maskmaxabund=1, numbands=None, band=None,
                        outfile=None, numthreads=1, logfile=sys.stderr):
    """
    Compute k-mer abundances for the specified sequence input.

    Expected input is a list of one or more FASTA/FASTQ files corresponding
    to a single sample. A counttable is created and populated with abundances
    of all k-mers observed in the input. If `mask` is provided, only k-mers not
    present in the mask will be loaded.
    """
    message = 'loading sample from ' + ','.join(seqfiles)
    print('[kevlar::counting]    ', message, file=logfile)

    sketch = khmer.Counttable(ksize, memory / 4, 4)
    n, nkmers = 0, 0
    for seqfile in seqfiles:
        parser = khmer.ReadParser(seqfile)
        threads = list()
        for _ in range(numthreads):
            if mask:
                if numbands:
                    thread = threading.Thread(
                        target=sketch.consume_seqfile_banding_with_mask,
                        args=(parser, numbands, band, mask, ),
                    )
                else:
                    thread = threading.Thread(
                        target=sketch.consume_seqfile_with_mask,
                        args=(parser, mask, ),
                    )
            else:
                if numbands:
                    thread = threading.Thread(
                        target=sketch.consume_seqfile_banding,
                        args=(parser, numbands, band, ),
                    )
                else:
                    thread = threading.Thread(
                        target=sketch.consume_seqfile,
                        args=(parser, ),
                    )
            threads.append(thread)
            thread.start()

    for thread in threads:
        thread.join()

    message = 'done loading reads'
    if numbands:
        message += ' (band {:d}/{:d})'.format(band+1, numbands)
    fpr = kevlar.sketch.estimate_fpr(sketch)
    message += '; {:d} reads processed'.format(parser.num_reads)
    message += ', {:d} distinct k-mers stored'.format(sketch.n_unique_kmers())
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    if fpr > maxfpr:
        message += ' (FPR too high, bailing out!!!)'
        raise SystemExit(message)

    if outfile:
        if not outfile.endswith(('.ct', '.counttable')):
            outfile += '.counttable'
        sketch.save(outfile)
        message += '; saved to "{:s}"'.format(outfile)
    print('[kevlar::counting]    ', message, file=logfile)

    return sketch


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
        sketch = load_sample_seqfile(
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
