#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict
import argparse
import re
import sys
import textwrap

import khmer
from khmer.utils import write_record
from khmer import khmer_args
import kevlar
import screed


class CaseSampleMismatchError(ValueError):
    pass


def subparser(subparsers):
    """Define the `kevlar novel` command-line interface."""

    desc = """\
    Identify "interesting" (potentially novel) k-mers and output the
    corresponding reads. Here we define "interesting" k-mers as those which are
    high abundance in each case sample and effectively absent (below some
    specified abundance threshold) in each control sample."""
    desc = textwrap.dedent(desc)
    epilog = """\
    Example::

        kevlar novel --out novel-reads.augfastq --case proband-reads.fq.gz \\
            --control father-reads-r1.fq.gz father-reads-r2.fq.gz \\
            --control mother-reads.fq.gz

    Example::

        kevlar novel --out novel-reads.augfastq.gz \\
            --control-counts father.counttable mother.counttable \\
            --case-counts proband.counttable --case proband-reads.fastq \\
            --ctrl-max 0 --case-min 10 --ksize 27

    Example::

        kevlar novel --out output.augfastq \\
            --case proband1.fq --case proband2.fq \\
            --control control1a.fq control1b.fq \\
            --control control2a.fq control2b.fq"""
    epilog = textwrap.dedent(epilog)
    subparser = subparsers.add_parser(
        'novel', description=desc, epilog=epilog, add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    samp_desc = """\
    Specify input files, as well as thresholds for selecting "interesting"
    k-mers. A single pass is made over input files for control samples (to
    compute k-mer abundances), while two passes are made over input files for
    case samples (to compute k-mer abundances, and then to identify
    "interesting" k-mers). The k-mer abundance computing steps can be skipped
    if pre-computed k-mer abunandances are provided using the "--case-counts"
    and/or "--control-counts" settings. If "--control-counts" is declared, then
    all "--control" flags are ignored. If "--case-counts" is declared,
    FASTA/FASTQ files must still be provided with "--case" for selecting
    "interesting" k-mers and reads."""
    samp_desc = textwrap.dedent(samp_desc)
    samp_args = subparser.add_argument_group('Case/control config', samp_desc)
    samp_args.add_argument(
        '--case', metavar='F', nargs='+', required=True, action='append',
        help='one or more FASTA/FASTQ files containing reads from a case '
        'sample; can be declared multiple times corresponding to multiple '
        'case samples, see examples below'
    )
    samp_args.add_argument(
        '--case-counts', metavar='F', nargs='+',
        help='counttable file(s) corresponding to each case sample; if not '
        'provided, k-mer abundances will be computed from FASTA/FASTQ input; '
        'only one counttable per sample, see examples below'
    )
    samp_args.add_argument(
        '--control', metavar='F', nargs='+', action='append',
        help='one or more FASTA/FASTQ files containing reads from a control '
        'sample; can be declared multiple times corresponding to multiple '
        'control samples, see examples below'
    )
    samp_args.add_argument(
        '--control-counts', metavar='F', nargs='+',
        help='counttable file(s) corresponding to each control sample; if not '
        'provided, k-mer abundances will be computed from FASTA/FASTQ input; '
        'only one counttable per sample, see examples below'
    )
    samp_args.add_argument(
        '-x', '--ctrl-max', metavar='X', type=int, default=1,
        help='k-mers with abund > X in any control sample are uninteresting; '
        'default is X=1'
    )
    samp_args.add_argument(
        '-y', '--case-min', metavar='Y', type=int, default=5,
        help='k-mers with abund < Y in any case sample are uninteresting; '
        'default is Y=5'
    )
    samp_args.add_argument(
        '-M', '--memory', default='1e6', type=khmer_args.memory_setting,
        metavar='MEM',
        help='total memory allocated to k-mer abundance for each sample; '
        'default is 1M; ignored when pre-computed k-mer abundances are '
        'supplied via counttable'
    )
    samp_args.add_argument(
        '--max-fpr', type=float, default=0.2, metavar='FPR',
        help='terminate if the expected false positive rate for any sample is '
        'higher than the specified FPR; default is 0.2'
    )

    band_desc = """\
    If memory is a limiting factor, it is possible to get a linear decrease in
    memory consumption by running `kevlar novel` in "banded" mode. Splitting
    the hashed k-mer space into N bands and only considering k-mers from one
    band at a time reduces the memory consumption to approximately 1/N of the
    total memory required. This implements a scatter/gather approach in which
    `kevlar novel` is run N times, after the results are combined using
    `kevlar filter`."""
    band_desc = textwrap.dedent(band_desc)
    band_args = subparser.add_argument_group('K-mer banding', band_desc)
    band_args.add_argument('--num-bands', type=int, metavar='N', default=0,
                           help='number of bands into which to divide the '
                           'hashed k-mer space')
    band_args.add_argument('--band', type=int, metavar='I', default=0,
                           help='a number between 1 and N (inclusive) '
                           'indicating the band to be processed')

    misc_args = subparser.add_argument_group('Miscellaneous settings')
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    misc_args.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    misc_args.add_argument('--upint', type=float, default=1e6, metavar='INT',
                           help='update interval for log messages; default is '
                           '1000000 (1 update message per millon reads)')


def load_sketch(infile, ksize, memory, max_fpr=0.2, numbands=0, band=0,
                logfile=sys.stderr):
    """Load a count table from disk."""
    message = 'Loading sample {}...'.format(infile)
    print('[kevlar::novel]    ', message, sep='', end='', file=logfile)

    if not ct.endswith(('.ct', '.counttable')):
        message = 'WARNING: counttable file does not have the expected '
        message += 'file extension; expect failures to occur'
        print('[kevlar::novel]        ', message, sep='', file=logfile)

    sketch = kevlar.sketch_autoload(
        infile, count=True, graph=False, ksize=ksize, table_size=memory/4,
        num_bands=numbands, band=band
    )

    fpr = kevlar.calc_fpr(sketch)
    message = 'done! estimated false positive rate is {:1.3f}'.format(fpr)
    if fpr > max_fpr:
        message += ' (FPR too high, bailing out!!!)'
        raise SystemExit(message)
    else:
        print(message, file=logfile)

    return sketch


def load_samples(samples, ksize, memory, counttables=None, max_fpr=0.2,
                 numbands=0, band=0, logfile=sys.stderr):
    sample_counts = list()
    if counttables:
        numcases = len(counttables)
        message = 'counttables for {:d} samples provided'.format(numcases)
        message += ', any corresponding FASTA/FASTQ input will be ignored '
        message += 'for computing k-mer abundances'
        print('[kevlar::novel] INFO:', message, file=logfile)
        for ct in counttables:
            table = load_sketch(ct, ksize, 1, max_fpr=max_fpr, logfile=logfile)
            sample_counts.append(table)
    else:
        for file_list in samples:
            counttable = kevlar.allocate_sketch(ksize, memory / 4, count=True)
            for fastx in file_list:
                if numbands > 1:
                    counttable.consume_seqfile_banding(fastx, numbands, band)
                else:
                    counttable.consume_seqfile(fastx)
            sample_counts.append(counttable)
    return sample_counts


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


def iter_read_multi_file(filenames):
    for filename in filenames:
        for record in screed.open(filename):
            yield record


def main(args):
    if (not args.num_bands) is not (not args.band):
        raise ValueError('Must specify --num-bands and --band together')

    timer = kevlar.Timer()
    timer.start()

    print('[kevlar::novel] Loading case samples', file=args.logfile)
    timer.start('loadall')
    timer.start('loadcases')
    cases = load_samples(
        args.case, args.ksize, args.memory, counttables=args.case_counts,
        max_fpr=args.max_fpr, numbands=args.num_bands, band=args.band-1,
        logfile=args.logfile
    )
    if args.case_counts:
        if len(cases) != args.case:
            message = '{:d} case samples declared '.format(len(args.cases))
            message += 'but {:d} counttables provided'.format(len(cases))
            raise CaseSampleMismatchError(message)
    elapsed = timer.stop('loadcases')
    print('[kevlar::novel] Case samples loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    timer.start('loadctrl')
    print('[kevlar::novel] Loading control samples', file=args.logfile)
    controls = load_samples(
        args.control, args.ksize, args.memory, counttables=args.control_counts,
        max_fpr=args.max_fpr, numbands=args.num_bands, band=args.band-1,
        logfile=args.logfile
    )
    elapsed = timer.stop('loadctrl')
    print('[kevlar::novel] Cntrl samples loaded in {:.2f} sec'.format(elapsed),
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
    for n, record in enumerate(iter_read_multi_file(infiles)):
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
