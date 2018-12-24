#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar


def filter(counts, reads, casemin=6, ctrlmax=1):
    for read in reads:
        validated_kmers = list()
        for ikmer in read.annotations:
            ikseq = read.ikmerseq(ikmer)
            ctrltoohigh = sum([1 for a in ikmer.abund[1:] if a > ctrlmax])
            if ctrltoohigh:
                continue

            newcount = counts.get(ikseq)
            casetoolow = newcount < minabund
            if casetoolow:
                continue
            ikmer.abund[0] = newcount
            validated_kmers.append(ikmer)
        if len(validated_kmers) == 0:
            continue
        yield read


def main(args):
    timer = kevlar.Timer()
    timer.start()

    counts = kevlar.sketch.load(args.counts)
    readstream = kevlar.parse_augmented_fastx(args.augfastq)
    outstream = kevlar.open(args.out, 'w')
    filterstream = filter(
        counts, readstream, casemin=args.case_min, ctrlmax=args.ctrl_max,
    )
    for n, record in enumerate(filterstream):
        kevlar.print_augmented_fastx(record, outstream)

    total = timer.stop()
    message = 'Total time: {:d} reads in {:.2f} seconds'.format(n + 1, total)
    kevlar.plog('[kevlar::filter]', message)
