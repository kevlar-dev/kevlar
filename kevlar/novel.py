#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import re
import sys

import khmer
from khmer import khmer_args
import kevlar


class KevlarCaseSampleMismatchError(ValueError):
    pass


def kmer_is_interesting(kmer, casecounts, controlcounts, case_min=5,
                        ctrl_max=1, screen_thresh=None):
    """
    Well, is it?

    Consult the k-mer's abundance in each sample to determine whether it is
    "interesting". It must be >= `case_min` in each of `casecounts`, and must
    be <= `ctrl_max` in each of `controlcounts`.

    Returns 4 values: 2 booleans, and 2 lists of integers
    - boolean indicating whether the k-mer is interesting
    - boolean indicating whether the entire read should be discarded
    - list of case sample abundances (empty if k-mer not interesting)
    - list of control sample abundances (empty if k-mer not interesting)
    """
    caseabunds = list()
    for ct in casecounts:
        abund = ct.get(kmer)
        if abund < case_min:
            discard = False
            if screen_thresh and abund < screen_thresh:
                discard = True
            return False, discard, [], []
        caseabunds.append(abund)

    ctrlabunds = list()
    for count in controlcounts:
        abund = count.get(kmer)
        if abund > ctrl_max:
            return False, False, [], []
        ctrlabunds.append(abund)

    return True, False, caseabunds, ctrlabunds


def main(args):
    if (not args.num_bands) is not (not args.band):
        raise ValueError('Must specify --num-bands and --band together')
    myband = args.band - 1 if args.band else None

    timer = kevlar.Timer()
    timer.start()

    timer.start('loadall')
    print('[kevlar::novel] Loading control samples', file=args.logfile)
    timer.start('loadctrl')
    if args.control_counts:
        numctrls = len(args.control_counts)
        message = 'counttables for {:d} sample(s) provided'.format(numctrls)
        message += ', any corresponding FASTA/FASTQ input will be ignored '
        message += 'for computing k-mer abundances'
        print('[kevlar::novel]    INFO:', message, file=args.logfile)
        controls = kevlar.counting.load_samples_sketchfiles(
            args.control_counts, args.max_fpr, args.logfile
        )
    else:
        controls = kevlar.counting.load_samples(
            args.control, args.ksize, args.memory, maxfpr=args.max_fpr,
            memfraction=None, numbands=args.num_bands, band=myband,
            logfile=args.logfile
        )
    elapsed = timer.stop('loadctrl')
    message = 'Control samples loaded in {:.2f} sec'.format(elapsed)
    print('[kevlar::novel]', message, file=args.logfile)

    print('[kevlar::novel] Loading case samples', file=args.logfile)
    timer.start('loadcases')
    if args.case_counts:
        numcases = len(args.case_counts)
        if numcases != len(args.case):
            message = '{:d} case samples declared '.format(len(args.case))
            message += 'but {:d} counttables provided'.format(numcases)
            raise KevlarCaseSampleMismatchError(message)
        message = 'counttables for {:d} samples provided'.format(numcases)
        message += ', any corresponding FASTA/FASTQ input will be ignored '
        message += 'for computing k-mer abundances'
        print('[kevlar::novel]    INFO:', message, file=args.logfile)
        cases = kevlar.counting.load_samples_sketchfiles(
            args.case_counts, args.max_fpr, args.logfile
        )
    else:
        cases = kevlar.counting.load_samples(
            args.case, args.ksize, args.memory, maxfpr=args.max_fpr,
            memfraction=None, numbands=args.num_bands, band=myband,
            logfile=args.logfile
        )
    elapsed = timer.stop('loadcases')
    print('[kevlar::novel] Case samples loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)
    elapsed = timer.stop('loadall')
    print('[kevlar::novel] All samples loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    timer.start('iter')
    ncases = len(args.case)
    message = 'Iterating over reads from {:d} case sample(s)'.format(ncases)
    print('[kevlar::novel]', message, file=args.logfile)
    nkmers = 0
    nreads = 0
    unique_kmers = set()
    outstream = kevlar.open(args.out, 'w')
    infiles = [f for filelist in args.case for f in filelist]
    for n, record in enumerate(kevlar.multi_file_iter_screed(infiles), 1):
        if n > 0 and n % args.upint == 0:
            elapsed = timer.probe('iter')
            msg = '    processed {} reads'.format(n)
            msg += ' in {:.2f} seconds...'.format(elapsed)
            print(msg, file=args.logfile)
        if len(record.sequence) < args.ksize:
            continue
        if re.search('[^ACGT]', record.sequence):
            # This check should be temporary; hopefully khmer will handle
            # this soon.
            continue

        discard_read = False
        record.ikmers = list()
        for i, kmer in enumerate(cases[0].get_kmers(record.sequence)):
            if args.num_bands:
                khash = cases[0].hash(kmer)
                if khash & (args.num_bands - 1) != args.band - 1:
                    continue
            interesting, discard, caseabund, ctrlabund = kmer_is_interesting(
                kmer, cases, controls, case_min=args.case_min,
                ctrl_max=args.ctrl_max, screen_thresh=args.abund_screen,
            )
            if discard:
                discard_read = True
                break
            if not interesting:
                continue
            abund = caseabund + ctrlabund
            ikmer = kevlar.KmerOfInterest(sequence=kmer, offset=i, abund=abund)
            record.ikmers.append(ikmer)
            minkmer = kevlar.revcommin(kmer)
            unique_kmers.add(minkmer)

        read_kmers = len(record.ikmers)
        if discard_read or read_kmers == 0:
            continue

        nreads += 1
        nkmers += read_kmers
        kevlar.print_augmented_fastx(record, outstream)

    elapsed = timer.stop('iter')
    message = 'Iterated over {} reads in {:.2f} seconds'.format(n, elapsed)
    print('[kevlar::novel]', message, file=args.logfile)

    message = 'Found {:d} instances'.format(nkmers)
    message += ' of {:d} unique novel kmers'.format(len(unique_kmers))
    message += ' in {:d} reads'.format(nreads)
    print('[kevlar::novel]', message, file=args.logfile)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::novel]', message, file=args.logfile)
