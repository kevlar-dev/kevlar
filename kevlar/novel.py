#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict
import argparse
import re
import sys

import khmer
from khmer.utils import write_record
from khmer import khmer_args
import kevlar
import screed


def subparser(subparsers):
    subparser = subparsers.add_parser('novel', add_help=False)

    samp_args = subparser.add_argument_group(
        'Case and control configuration',
        'Specify input files and thresholds for identifying "interesting" '
        'k-mers, which are high abundance in each case sample and absent (or '
        'low abundance) in each control sample.'
    )
    samp_args.add_argument('--cases', metavar='F', nargs='+',
                           required=True, help='one or more Fastq files '
                           'corresponding to case sample(s)')
    samp_args.add_argument('--controls', metavar='F', nargs='+',
                           required=True, help='one or more Fastq files '
                           'corresponding to control sample(s)')
    samp_args.add_argument('-x', '--ctrl_max', metavar='X', type=int,
                           default=1, help='k-mers with abund > X in any '
                           'control sample are uninteresting; default=1')
    samp_args.add_argument('-y', '--case_min', metavar='Y',
                           type=int, default=5, help='ignore k-mers from case '
                           'with abund < Y; default=5')
    samp_args.add_argument('-M', '--memory', default='1e6',
                           type=khmer_args.memory_setting, metavar='MEM',
                           help='total memory to allocate for each sample; '
                           'default is 1M')
    samp_args.add_argument('--max-fpr', type=float, default=0.2, metavar='FPR',
                           help='terminate if the expected false positive rate'
                           ' for any sample is higher than the specified FPR; '
                           'default is 0.2')

    band_args = subparser.add_argument_group(
        'K-mer banding',
        'If memory is a limiting factor, it is possible to get a linear '
        'decrease in memory consumption by running `kevlar novel` in "banded" '
        'mode. Splitting the hashed k-mer space into N bands and only '
        'considering k-mers from one band at a time reduces the memory '
        'consumption to approximately 1/N of the total memory required. This '
        'implements a scatter/gather approach in which `kevlar novel` is run N'
        ' times, after the results are combined using `kevlar filter`.'
    )
    band_args.add_argument('--num-bands', type=int, metavar='N', default=None,
                           help='number of bands into which to divide the '
                           'hashed k-mer space')
    band_args.add_argument('--band', type=int, metavar='I', default=None,
                           help='a number between 1 and N (inclusive) '
                           'indicating the band to be processed')

    misc_args = subparser.add_argument_group(
        'Miscellaneous settings'
    )
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    misc_args.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    misc_args.add_argument('--upint', type=float, default=1e6, metavar='INT',
                           help='update interval for log messages; default is '
                           '1000000 (1 update message per millon reads)')


def load_samples(samples, ksize, memory, max_fpr=0.2, numbands=None, band=None,
                 logfile=sys.stderr):
    tables = list()
    for sample in samples:
        print('[kevlar::novel]     Loading counttable', sample, '...', end='',
              file=logfile)
        ct = khmer.Counttable(ksize, memory / 4, 4)
        if numbands:
            nr, nk = ct.consume_seqfile_banding(sample, numbands, band - 1)
            message = 'batch {:d}/{:d} done'.format(band, numbands)
        else:
            message = 'done'
            nr, nk = ct.consume_seqfile(sample)
        fpr = kevlar.calc_fpr(ct)
        message += ', k={:d}'.format(ct.ksize())
        message += '; {:d} reads and {:d} k-mers consumed'.format(nr, nk)
        message += '; estimated false positive rate is {:1.3f}'.format(fpr)
        if fpr > max_fpr:
            raise SystemExit(message)
        else:
            print(message, file=logfile)
        tables.append(ct)

    return tables


def kmer_is_interesting(kmer, casecounts, controlcounts, case_min=5,
                        ctrl_max=1):
    caseabunds = list()
    for ct in casecounts:
        abund = ct.get(kmer)
        if abund < case_min:
            return False
        caseabunds.append(abund)

    ctrlabunds = list()
    for count in controlcounts:
        abund = count.get(kmer)
        if abund > ctrl_max:
            return False
        ctrlabunds.append(abund)

    return caseabunds + ctrlabunds


def iter_screed(filenames):
    for filename in filenames:
        for record in screed.open(filename):
            yield record


def main(args):
    if (args.num_bands is None) is not (args.band is None):
        raise ValueError('Must specify --num-bands and --band together')

    timer = kevlar.Timer()
    timer.start()

    print('[kevlar::novel] Loading case samples', file=args.logfile)
    timer.start('loadall')
    timer.start('loadcase')
    cases = load_samples(args.cases, args.ksize, args.memory, args.max_fpr,
                         args.num_bands, args.band, args.logfile)
    elapsed = timer.stop('loadcase')
    print('[kevlar::novel] Case samples loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    elapsed = timer.start('loadctrl')
    print('[kevlar::novel] Loading control samples', file=args.logfile)
    controls = load_samples(args.controls, args.ksize, args.memory,
                            args.max_fpr, args.num_bands, args.band,
                            args.logfile)
    elapsed = timer.stop('loadctrl')
    print('[kevlar::novel] Cntrl samples loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)
    elapsed = timer.stop('loadall')
    print('[kevlar::novel] All samples loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    timer.start('iter')
    print('[kevlar::novel] Iterating over case reads', args.cases,
          file=args.logfile)
    nkmers = 0
    nreads = 0
    unique_kmers = set()
    outstream = kevlar.open(args.out, 'w')
    for n, record in enumerate(iter_screed(args.cases)):
        if n > 0 and n % args.upint == 0:
            elapsed = timer.probe('iter')
            msg = '    processed {} reads'.format(n)
            msg += ' in {:.2f} seconds...'.format(elapsed)
            print(msg, file=args.logfile)
        if re.search('[^ACGT]', record.sequence):
            # This check should be temporary; hopefully khmer will handle
            # this soon.
            continue

        record.ikmers = list()
        for i, kmer in enumerate(cases[0].get_kmers(record.sequence)):
            if args.num_bands:
                khash = cases[0].hash(kmer)
                if khash & (args.num_bands - 1) != args.band - 1:
                    continue
            abund = kmer_is_interesting(kmer, cases, controls, args.case_min,
                                        args.ctrl_max)
            if not abund:
                continue
            ikmer = kevlar.KmerOfInterest(sequence=kmer, offset=i, abund=abund)
            record.ikmers.append(ikmer)
            minkmer = kevlar.revcommin(kmer)
            unique_kmers.add(minkmer)

        read_kmers = len(record.ikmers)
        if read_kmers > 0:
            nreads += 1
            nkmers += read_kmers
            kevlar.print_augmented_fastq(record, outstream)

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
