#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
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


def subparser(subparsers):
    subparser = subparsers.add_parser('count', add_help=False)

    samp_args = subparser.add_argument_group(
        'Case and control configuration',
        'Specify input files for k-mer counting. The stored k-mer counts will '
        'subsequently be used to identify "interesting k-mers", or those '
        'k-mers that are high abundance in the case and absent (or at least '
        'low abundance) in controls. Here we specify the abundance threshold '
        'at which a k-mer can be considered effectively "absent" from a '
        'control sample.'
    )
    samp_args.add_argument('--case', metavar='F', required=True,
                           help='one or more Fastq files corresponding to '
                           'case sample(s)')
    samp_args.add_argument('--controls', metavar='F', nargs='+',
                           required=True, help='one or more Fastq files '
                           'corresponding to control sample(s)')
    samp_args.add_argument('-x', '--ctrl_max', metavar='X', type=int,
                           default=1, help='k-mers with abund > X in any '
                           'control sample are uninteresting; default=1')

    memory_args = subparser.add_argument_group(
        'Memory allocation',
        'Specify how much memory to allocate for the sketch data structures '
        'used to store k-mer counts. The first control sample will be '
        'allocated the full amount of specifed `--memory`, while the case '
        'sample and all subsequent control samples will be allocated a '
        'fraction thereof.'
    )
    memory_args.add_argument('-M', '--memory', default='1e6',
                             type=khmer_args.memory_setting, metavar='MEM',
                             help='total memory to allocate for the initial '
                             'control sample; default is 1M')
    memory_args.add_argument('-f', '--mem-frac', type=float, default=0.1,
                             metavar='F', help='fraction of the total memory '
                             'to allocate to subsequent samples; default is '
                             '0.1')
    memory_args.add_argument('--max-fpr', type=float, default=0.2,
                             metavar='FPR', help='terminate if the expected '
                             'false positive rate for any sample is higher '
                             'than the specified FPR; default is 0.2')

    band_args = subparser.add_argument_group(
        'K-mer banding',
        'If memory is a limiting factor, it is possible to get a linear '
        'decrease in memory consumption by running `kevlar novel` in "banded" '
        'mode. Splitting the hashed k-mer space into N bands and only '
        'considering k-mers from one band at a time reduces the memory '
        'consumption to approximately 1/N of the total memory required. This '
        'implements a scatter/gather approach in which `kevlar novel` is run N'
        ' times, after which the results are combined using `kevlar filter`.'
    )
    band_args.add_argument('--num-bands', type=int, metavar='N', default=0,
                           help='number of bands into which to divide the '
                           'hashed k-mer space')
    band_args.add_argument('--band', type=int, metavar='I', default=0,
                           help='a number between 1 and N (inclusive) '
                           'indicating the band to be processed')

    misc_args = subparser.add_argument_group(
        'Miscellaneous settings'
    )
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')


def load_controls(samples, ksize, memory, memfraction=0.1, max_fpr=0.2,
                  maxabund=1, numbands=None, band=None, logfile=sys.stderr):
    tables = list()
    for n, sample in enumerate(samples):
        message = 'loading control {:d} from "{:s}"'.format(n + 1, sample)
        print('[kevlar::count]    ', message, file=logfile)
        sketchmem = memory if n == 0 else memory * memfraction
        ct = khmer.Counttable(ksize, sketchmem / 4, 4)
        reads = khmer.ReadParser(sample)
        kmers_stored = 0
        for i, record in enumerate(reads):
            for j, kmer in enumerate(ct.get_kmers(record.sequence)):
                if re.search('[^ACGT]', kmer):
                    continue
                if numbands:
                    khash = ct.hash(kmer)
                    if khash & (numbands - 1) != band - 1:
                        continue
                if n == 0:
                    ct.add(kmer)
                    kmers_stored += 1
                else:
                    for table in tables:
                        if table.get(kmer) > maxabund:
                            break
                    else:
                            ct.add(kmer)
                            kmers_stored += 1

        message = 'done loading control {:d}'.format(n + 1)
        if numbands:
            message += ' (band {:d}/{:d})'.format(band, numbands)
        fpr = kevlar.calc_fpr(ct)
        message += '; {:d} reads processed'.format(i + 1)
        message += ', {:d} k-mers stored'.format(kmers_stored)
        message += '; estimated false positive rate is {:1.3f}'.format(fpr)
        if fpr > max_fpr:
            message += ' (FPR too high, bailing out!!!)'
            raise SystemExit(message)
        else:
            print('[kevlar::count]    ', message, file=logfile)

        tables.append(ct)
    return tables


def load_case(sample, controls, ksize, memory, memfraction=0.1, max_fpr=0.2,
              ctrl_maxabund=1, numbands=None, band=None, logfile=sys.stderr):
    message = 'loading case from "{:s}"'.format(sample)
    print('[kevlar::count]    ', message, file=logfile)
    sketchmem = memory * memfraction
    ct = khmer.Counttable(ksize, sketchmem / 4, 4)
    reads = khmer.ReadParser(sample)
    kmers_stored = 0
    for i, record in enumerate(reads):
        for j, kmer in enumerate(ct.get_kmers(record.sequence)):
            if re.search('[^ACGT]', kmer):
                continue
            if numbands:
                khash = ct.hash(kmer)
                if khash & (numbands - 1) != band - 1:
                    continue
            for ctrl in controls:
                if ctrl.get(kmer) > ctrl_maxabund:
                    break
            else:
                ct.add(kmer)
                kmers_stored += 1

    message = 'done loading case sample'
    if numbands:
        message += ' (band {:d}/{:d})'.format(band, numbands)
    fpr = kevlar.calc_fpr(ct)
    message += '; {:d} reads processed'.format(i + 1)
    message += ', {:d} k-mers stored'.format(kmers_stored)
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    if fpr > max_fpr:
        message += ' (FPR too high, bailing out!!!)'
        raise SystemExit(message)
    else:
        print('[kevlar::count]    ', message, file=logfile)

    return ct


def table_filename(sample, band=0):
    suffixes = ['.fa', '.fasta', '.fna', '.fq', '.fastq']
    gzipsuffixes = [s + '.gz' for s in suffixes]
    if sample.endswith(tuple(suffixes)):
        prefix = '.'.join(sample.split('.')[:-1])
    elif sample.endswith(tuple(gzipsuffixes)):
        prefix = '.'.join(sample.split('.')[:-2])
    else:
        prefix = sample
    if band:
        prefix += '.band{:d}'.format(band)
    return prefix + '.counttable'


def save_tables(case, casetable, controls, controltables, band=0):
    casename = table_filename(case, band)
    casetable.save(casename)

    for control, controltable in zip(controls, controltables):
        controlname = table_filename(control)
        controltable.save(controlname)


def main(args):
    if (args.num_bands is None) is not (args.band is None):
        raise ValueError('Must specify --num-bands and --band together')

    timer = kevlar.Timer()
    timer.start()

    timer.start('loadall')
    timer.start('loadctrl')
    print('[kevlar::count] Loading control samples', file=args.logfile)
    controls = load_controls(args.controls, args.ksize, args.memory,
                             args.mem_frac, args.max_fpr, args.ctrl_max,
                             args.num_bands, args.band, args.logfile)
    elapsed = timer.stop('loadctrl')
    print('[kevlar::count] Cntrl samples loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)
    elapsed = timer.stop('loadall')
    print('[kevlar::count] All samples loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    print('[kevlar::count] Loading case sample', file=args.logfile)
    timer.start('loadcase')
    case = load_case(args.case, controls, args.ksize, args.memory,
                     args.mem_frac, args.max_fpr, args.ctrl_max,
                     args.num_bands, args.band, args.logfile)
    elapsed = timer.stop('loadcase')
    print('[kevlar::count] Case samples loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    print('[kevlar::count] Saving count tables', file=args.logfile)
    timer.start('savetables')
    save_tables(args.case, case, args.controls, controls, args.band)
    elapsed = timer.stop('savetables')
    print('[kevlar::count] Count tables saved in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::count]', message, file=args.logfile)
