#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

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


def load_samples(counttables=None, filelists=None, ksize=31, memory=1e6,
                 maxfpr=0.2, numbands=None, band=None, numthreads=1,
                 logstream=sys.stderr):
    assert counttables or filelists
    if counttables:
        numctrls = len(counttables)
        message = 'counttables for {:d} sample(s) provided'.format(numctrls)
        message += ', any corresponding FASTA/FASTQ input will be ignored '
        message += 'for computing k-mer abundances'
        print('[kevlar::novel]    INFO:', message, file=logstream)
        samples = kevlar.sketch.load_sketchfiles(
            counttables, maxfpr, logstream,
        )
    else:
        samples = list()
        for filelist in filelists:
            sample = kevlar.count.load_sample_seqfile(
                filelist, ksize, memory, maxfpr=maxfpr, numbands=numbands,
                band=band, numthreads=numthreads, logfile=logstream
            )
            samples.append(sample)
    return samples


def save_counts(filelist, tablelist, logstream=sys.stderr):
    if len(filelist) != len(tablelist):
        msg = 'number of filenames provided ({:d})'.format(len(filelist))
        msg += 'does not match the number of '
        msg += 'samples provided ({:d})'.format(len(tablelist))
        msg += '; stubbornly refusing to save k-mer counts'
        print('[kevlar::novel] WARNING:', msg, file=logstream)
        return
    for outfile, counttable in zip(filelist, tablelist):
        if not outfile.endswith(('.ct', '.counttable')):
            outfile += '.counttable'
        print('   ', outfile, file=logstream)
        counttable.save(outfile)


def save_all_counts(casecounts, casefiles, ctrlcounts, ctrlfiles,
                    logstream=sys.stderr):
    if casefiles:
        message = 'saving k-mer counts for case samples'
        print('[kevlar::novel]', message, file=logstream)
        save_counts(casefiles, casecounts, logstream=logstream)
    if ctrlfiles:
        message = 'saving k-mer counts for control samples'
        print('[kevlar::novel]', message, file=logstream)
        save_counts(ctrlfiles, ctrlcounts, logstream=logstream)


def novel(casestream, casecounts, controlcounts, ksize=31, abundscreen=None,
          casemin=5, ctrlmax=0, numbands=None, band=None, skipuntil=None,
          logstream=sys.stderr):
    numbands_unset = not numbands
    band_unset = not band and band != 0
    if numbands_unset is not band_unset:
        raise ValueError('Must specify `numbands` and `band` together')

    if band is not None and band < 0:
        maxband = numbands - 1
        message = '`band` must be a value between 0 and {:d}'.format(maxband)
        message += ' (`numbands` - 1), inclusive'
        raise ValueError(message)

    timer = kevlar.Timer()
    timer.start()
    nkmers = 0
    nreads = 0
    update_message = '[kevlar::novel]     processed {counter} reads'
    skip_message = None
    if skipuntil:
        msg = '; skipping reads in search of {read}'.format(read=skipuntil)
        skip_message = update_message + msg
    first_message = skip_message if skipuntil else update_message
    progress_indicator = kevlar.ProgressIndicator(
        first_message, interval=1e6,
        breaks=[1e7, 1e8, 1e9], usetimer=True, logstream=logstream,
    )
    unique_kmers = set()
    for n, record, mate in kevlar.paired_reader(casestream):
        progress_indicator.update()
        if skipuntil:  # pragma: no cover
            if record.name == skipuntil:
                message = 'Found read {:s}'.format(skipuntil)
                message += ' (skipped {:d} reads)'.format(n)
                print('[kevlar::novel]', message, file=logstream)
                skipuntil = False
                progress_indicator.message = update_message
            continue

        if len(record.sequence) < ksize:
            continue
        if re.search('[^ACGT]', record.sequence):
            # This check should be temporary; hopefully khmer will handle
            # this soon.
            continue

        discard_read = False
        irecord = None
        for i, kmer in enumerate(casecounts[0].get_kmers(record.sequence)):
            if numbands:
                khash = casecounts[0].hash(kmer)
                if khash & (numbands - 1) != band - 1:
                    continue
            interesting, discard, caseabund, ctrlabund = kmer_is_interesting(
                kmer, casecounts, controlcounts, case_min=casemin,
                ctrl_max=ctrlmax, screen_thresh=abundscreen,
            )
            if discard:
                discard_read = True
                break
            if not interesting:
                continue
            if irecord is None:
                irecord = kevlar.sequence.copy_record(record)
            abund = tuple(caseabund + ctrlabund)
            irecord.annotate(kmer, i, abund)
            minkmer = kevlar.revcommin(kmer)
            unique_kmers.add(minkmer)

        if discard_read or irecord is None:
            continue

        nreads += 1
        nkmers += len(irecord.annotations)
        if mate:
            irecord.add_mate(mate.sequence)
        yield irecord

    elapsed = timer.stop()
    message = 'Found {:d} instances'.format(nkmers)
    message += ' of {:d} unique novel kmers'.format(len(unique_kmers))
    message += ' in {:d} reads'.format(nreads)
    message += ' in {:.2f} seconds'.format(elapsed)
    print('[kevlar::novel]', message, file=logstream)


def main(args):
    timer = kevlar.Timer()
    timer.start()
    if (not args.num_bands) is not (not args.band):
        raise ValueError('Must specify --num-bands and --band together')
    myband = args.band - 1 if args.band else None

    timer.start('loadall')
    print('[kevlar::novel] Loading control samples', file=args.logfile)
    timer.start('loadctrl')
    controls = load_samples(
        args.control_counts, args.control, args.ksize, args.memory,
        args.max_fpr, args.num_bands, myband, args.threads, args.logfile
    )
    elapsed = timer.stop('loadctrl')
    message = 'Control samples loaded in {:.2f} sec'.format(elapsed)
    print('[kevlar::novel]', message, file=args.logfile)

    print('[kevlar::novel] Loading case samples', file=args.logfile)
    timer.start('loadcases')
    cases = load_samples(
        args.case_counts, args.case, args.ksize, args.memory,
        args.max_fpr, args.num_bands, myband, args.threads, args.logfile
    )
    elapsed = timer.stop('loadcases')
    print('[kevlar::novel] Case samples loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)
    elapsed = timer.stop('loadall')
    print('[kevlar::novel] All samples loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    save_all_counts(cases, args.save_case_counts, controls,
                    args.save_ctrl_counts, logstream=args.logfile)

    timer.start('iter')
    ncases = len(args.case)
    message = 'Iterating over reads from {:d} case sample(s)'.format(ncases)
    print('[kevlar::novel]', message, file=args.logfile)
    outstream = kevlar.open(args.out, 'w')
    infiles = [f for filelist in args.case for f in filelist]
    caserecords = kevlar.multi_file_iter_screed(infiles)
    readstream = novel(
        caserecords, cases, controls, ksize=args.ksize,
        abundscreen=args.abund_screen, casemin=args.case_min,
        ctrlmax=args.ctrl_max, numbands=args.num_bands, band=myband,
        skipuntil=args.skip_until, logstream=args.logfile,
    )
    for augmented_read in readstream:
        kevlar.print_augmented_fastx(augmented_read, outstream)

    elapsed = timer.stop('iter')
    message = 'Iterated over all case reads in {:.2f} seconds'.format(elapsed)
    print('[kevlar::novel]', message, file=args.logfile)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::novel]', message, file=args.logfile)
