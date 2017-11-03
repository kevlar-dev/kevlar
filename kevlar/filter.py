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

import khmer
from khmer import khmer_args
import kevlar


def load_mask(maskfiles, ksize, memory, maxfpr=0.001, savefile=None,
              logstream=sys.stderr):
    """Load reference genome and/or contaminant database from a file."""
    if maskfiles is None:
        return None

    timer = kevlar.Timer()
    timer.start('loadmask')
    print('[kevlar::filter] Loading mask from', maskfiles, file=logstream)

    if len(maskfiles) == 1 and maskfiles[0].endswith(('.nt', '.nodetable')):
        mask = kevlar.sketch.load(maskfiles[0])
        message = '    nodetable loaded'
    else:
        buckets = memory * khmer._buckets_per_byte['nodegraph'] / 4
        mask = khmer.Nodetable(ksize, buckets, 4)
        nr, nk = 0, 0
        for maskfile in maskfiles:
            numreads, numkmers = mask.consume_seqfile(maskfile)
            nr += numreads
            nk += numkmers
        message = '    {:d} sequences and {:d} k-mers consumed'.format(nr, nk)
    fpr = kevlar.sketch.estimate_fpr(mask)
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    print(message, file=logstream)
    if fpr > maxfpr:
        raise SystemExit('FPR too high, bailing out')
    if savefile:
        mask.save(savefile)
        message = '    nodetable saved to "{:s}"'.format(savefile)
        print(message, file=logstream)

    elapsed = timer.stop('loadmask')
    print('[kevlar::filter]', 'Mask loaded in {:.2f} sec'.format(elapsed),
          file=logstream)
    return mask


def summarize_readset(readset, logfile):
    fpr = kevlar.sketch.estimate_fpr(readset._counts)
    message = '    {:d} instances'.format(readset.read_instances)
    message += ' of {:d} reads consumed,\n'.format(readset.distinct_reads)
    message += '    annotated with'
    message += ' {:d} instances '.format(readset.ikmer_instances)
    message += 'of {:d} distinct'.format(readset.distinct_ikmers)
    message += ' "interesting" k-mers;\n'
    message += '    estimated false positive rate is {:1.3f}'.format(fpr)
    if logfile is not None:
        print(message, file=logfile)
    return fpr


def summarize_validate(readset, n, logfile=sys.stderr):
    int_distinct = readset.masked[0] + readset.lowabund[0] + readset.valid[0]
    int_instances = readset.masked[1] + readset.lowabund[1] + readset.valid[1]

    message = '    processed {:d} instances'.format(int_instances)
    message += ' of {:d} distinct "interesting" k-mers'.format(int_distinct)
    message += ' in {:d} reads'.format(len(readset))
    message += '\n        '
    message += '{:d} instances'.format(readset.masked[1])
    message += ' of {:d} distinct k-mers'.format(readset.masked[0])
    message += ' masked by the reference genome'
    message += '\n        '
    message += '{:d} instances'.format(readset.lowabund[1])
    message += ' of {:d} distinct k-mers'.format(readset.lowabund[0])
    message += ' discarded due to low abundance'
    message += '\n        '
    message += '{:d} instances'.format(readset.valid[1])
    message += ' of {:d} distinct k-mers'.format(readset.valid[0])
    message += ' validated as novel'
    message += '\n        '
    message += '{:d} reads'.format(readset.discarded)
    message += ' with no surviving valid k-mers ignored'
    message += '\n        '
    message += '{:d} reads written to output'.format(n)
    if logfile is not None:
        print(message, file=logfile)


def filter(readstream, mask=None, minabund=5, ksize=31, memory=1e6,
           maxfpr=0.001, logstream=sys.stderr):
    timer = kevlar.Timer()
    timer.start('recalc')
    print('[kevlar::filter] Loading input; recalculate k-mer abundances,',
          'de-duplicate reads and merge k-mers',
          file=logstream)
    readset = kevlar.seqio.AnnotatedReadSet(ksize, memory)
    for record in readstream:
        readset.add(record)
    fpr = summarize_readset(readset, logstream)
    if fpr > maxfpr:
        print('[kevlar::filter] FPR too high, bailing out', file=logstream)
        sys.exit(1)
    elapsed = timer.stop('recalc')
    print('[kevlar::filter] Input loaded in {:.2f} sec'.format(elapsed),
          file=logstream)

    timer.start('validate')
    print('[kevlar::filter] Validate k-mers and print reads',
          file=logstream)
    readset.validate(mask=mask, minabund=minabund)
    for n, record in enumerate(readset, 1):
        yield record
    summarize_validate(readset, n, logstream)
    elapsed = timer.stop('validate')
    print('[kevlar::filter] k-mers validated and reads printed',
          'in {:.2f} sec'.format(elapsed), file=logstream)


def main(args):
    timer = kevlar.Timer()
    timer.start()

    mask = load_mask(
        args.mask, args.ksize, args.mask_memory, maxfpr=args.mask_max_fpr,
        savefile=args.save_mask, logstream=args.logfile
    )
    readstream = kevlar.seqio.afxstream(args.augfastq)
    outstream = kevlar.open(args.out, 'w')
    filterstream = filter(
        readstream, mask, minabund=args.min_abund, ksize=args.ksize,
        memory=args.abund_memory, maxfpr=args.abund_max_fpr,
        logstream=args.logfile
    )
    for record in filterstream:
        kevlar.print_augmented_fastx(record, outstream)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::filter]', message, file=args.logfile)
