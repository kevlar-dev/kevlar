#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
from kevlar.sequence import KmerOfInterest
import khmer


def first_pass(reads, mask, memory, timer):
    kevlar.plog('[kevlar::filter] First pass: re-counting k-mers')
    timer.start('firstpass')
    counts = None
    progress_indicator = kevlar.ProgressIndicator(
        '[kevlar::filter]     processed {counter} reads',
        interval=1e5, breaks=[1e6, 1e7],
    )
    for n, read in enumerate(reads, 1):
        progress_indicator.update()
        if len(read.annotations) == 0:
            continue
        if counts is None:
            ksize = read.annotations[0].ksize
            counts = khmer.Counttable(ksize, memory / 4, 4)
        for ikmer in read.annotations:
            ikseq = read.ikmerseq(ikmer)
            if mask and mask.get(ikseq) > 0:
                continue
            counts.add(ikseq)
    elapsed = timer.stop('firstpass')
    message = 'First pass complete!'
    message += ' Processed {:d} reads in {:.2f} seconds!'.format(n, elapsed)
    kevlar.plog('[kevlar::filter]', message)
    return counts


def check_fpr(counts, maxfpr):
    fpr = kevlar.sketch.estimate_fpr(counts)
    message = 'FPR for re-computed k-mer counts: {:1.3f}'.format(fpr)
    kevlar.plog('[kevlar::filter]', message)
    if fpr > maxfpr:
        message += 'FPR too high, bailing out!!!'
        raise kevlar.sketch.KevlarUnsuitableFPRError(message)


def second_pass(reads, counts, casemin, ctrlmax, timer):
    kevlar.plog('[kevlar::filter] Second pass: discarding k-mers/reads')
    timer.start('secondpass')
    kept = 0
    progress_indicator = kevlar.ProgressIndicator(
        '[kevlar::filter]     processed {counter} reads',
        interval=1e5, breaks=[1e6, 1e7],
    )
    for read in reads:
        progress_indicator.update()
        validated_kmers = list()
        for ikmer in read.annotations:
            ikseq = read.ikmerseq(ikmer)
            ctrltoohigh = sum([1 for a in ikmer.abund[1:] if a > ctrlmax])
            if ctrltoohigh:
                continue
            newcount = counts.get(ikseq)
            casetoolow = newcount < casemin
            if casetoolow:
                continue
            newabund = tuple([newcount] + list(ikmer.abund[1:]))
            newikmer = KmerOfInterest(ikmer.ksize, ikmer.offset, newabund)
            validated_kmers.append(newikmer)
        if len(validated_kmers) == 0:
            continue
        yield read
        kept += 1
    elapsed = timer.stop('secondpass')
    message = 'Second pass complete!'
    message += ' Validated {:d} reads in {:.2f} seconds!'.format(kept, elapsed)
    kevlar.plog('[kevlar::filter]', message)


def filter(readfile, mask=None, memory=1e6, maxfpr=0.01, casemin=6, ctrlmax=1):
    timer = kevlar.Timer()
    timer.start()
    reader = kevlar.parse_augmented_fastx(kevlar.open(readfile, 'r'))
    counts = first_pass(reader, mask, memory, timer)
    check_fpr(counts, maxfpr)
    reader = kevlar.parse_augmented_fastx(kevlar.open(readfile, 'r'))
    for read in second_pass(reader, counts, casemin, ctrlmax, timer):
        yield read
    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    kevlar.plog('[kevlar::filter]', message)


def main(args):
    mask = kevlar.sketch.load(args.mask)
    with kevlar.open(args.out, 'w') as outstream:
        filterstream = filter(
            args.augfastq, mask=mask, memory=args.memory, maxfpr=args.max_fpr,
            casemin=args.case_min, ctrlmax=args.ctrl_max,
        )
        for n, record in enumerate(filterstream):
            kevlar.print_augmented_fastx(record, outstream)
