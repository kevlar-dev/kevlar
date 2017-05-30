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
    samp_args.add_argument('case', nargs='+',
                           help='one or more Fastq files containing reads from'
                           ' the case sample; if no counttables are provided '
                           'using the "--case-tab" setting, then "kevlar '
                           'novel" makes two passes over the files, one to '
                           'compute k-mer abundances, and a second to identify'
                           ' novel k-mers; if counttables are specified, the '
                           'Fastq files are used only to identify novel k-mers'
                           ' in a single pass')
    samp_args.add_argument('--case-counts', metavar='F',
                           help='counttable file corresponding to the case '
                           'sample; if not provided, will be computed from '
                           'Fastq input')
    samp_args.add_argument('--controls', metavar='F', nargs='+', required=True,
                           help='one or more Fastq files, each corresponding '
                           'to a distinct control sample; alternatively, '
                           'counttable files with k-mer abundances '
                           'pre-computed by "kevlar count" may be provided')
    samp_args.add_argument('-x', '--ctrl-max', metavar='X', type=int,
                           default=1, help='k-mers with abund > X in any '
                           'control sample are uninteresting; default=1')
    samp_args.add_argument('-y', '--case-min', metavar='Y',
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


def load_sample(sample, ksize, memory, max_fpr=0.2, numbands=0, band=0,
                logfile=sys.stderr):
    print('[kevlar::novel]     Loading sample ', sample, '...', sep='', end='',
          file=logfile)
    sketch = kevlar.sketch_autoload(
        sample, count=True, graph=False, ksize=ksize, tablesize=memory/4,
        num_bands=numbands, band=band
    )
    fpr = kevlar.calc_fpr(ct)
    message = 'done! estimated false positive rate is {:1.3f}'.format(fpr)
    if fpr > max_fpr:
        message += ' (FPR too high, bailing out!!!)'
        raise SystemExit(message)
    else:
        print(message, file=logfile)
    return sketch


def load_case(fastqs, ksize, memory, ct=None, max_fpr=0.2, numbands=0, band=0,
              logfile=sys.stderr):
    if ct:
        if not args.case_counts.endswith(('.ct', '.counttable')):
            message = 'counttable does not have the expected file extension;'
            message += ' expect failures to occur'
            print('[kevlar::novel] WARNING:', message, file=args.logfile)
        return load_sample(ct, ksize, 1, max_fpr=max_fpr, logfile=logfile)
    else:
        sketch = allocate_sketch(ksize, memory / 4, count=True)
        for fastq in fastqs:
            if numbands > 1:
                sketch.consume_seqfile_banding(fastq, numbands, band)
            else:
                sketch.consume_seqfile(fastq)
        return sketch


def load_controls(samples, ksize, memory, max_fpr=0.2, numbands=0, band=0,
                  logfile=sys.stderr):
    sketches = list()
    for sample in samples:
        sketch = load_sample(sample, ksize, memory, max_fpr=max_fpr,
                             numbands=numbands, band=band, logfile=logfile)
        sketches.append(sketch)
    return sketches


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
    case = load_case(
        args.case, args.ksize, args.memory, ct=args.case_counts,
        max_fpr=args.max_fpx, numbands=args.num_bands, band=args.band,
        logfile=sys.stderr
    )
    elapsed = timer.stop('loadcase')
    print('[kevlar::novel] Case samples loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    elapsed = timer.start('loadctrl')
    print('[kevlar::novel] Loading control samples', file=args.logfile)
    controls = load_controls(args.controls, args.ksize, args.memory,
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
    for n, record in enumerate(iter_screed(args.case)):
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
