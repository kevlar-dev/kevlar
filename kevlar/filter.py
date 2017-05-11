#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import argparse
import re
import sys

import khmer
from khmer import khmer_args
import kevlar


def subparser(subparsers):
    subparser = subparsers.add_parser('filter', add_help=False)
    subparser._positionals.title = 'Required inputs'

    refr_args = subparser.add_argument_group(
        'Reference genome',
        'A reference genome is not required, but if available can be supplied '
        'to discard "interesting" k-mers that turn out not to be novel.'
    )
    refr_args.add_argument('--refr', metavar='FILE', type=str, default=None,
                           help='reference genome in Fasta/Fastq format')
    refr_args.add_argument('--refr-memory', metavar='MEM', default='1e6',
                           type=khmer_args.memory_setting,
                           help='memory to allocate for storing the reference '
                           'genome; default is 1M')
    refr_args.add_argument('--refr-max-fpr', type=float, metavar='FPR',
                           default=0.001, help='terminate if the expected '
                           'false positive rate is higher than the specified '
                           'FPR; default is 0.001')

    contam_args = subparser.add_argument_group(
        'Screening for known contaminants',
        'It may not be possible to anticipate every possible contaminant that '
        'may be present in a sample, but if there is a common or known set of '
        'contaminants these can be filtered out at this step. Unknown '
        'contaminants will have to be filtered out after read partitioning.'
    )
    contam_args.add_argument('--contam', metavar='FILE', type=str,
                             default=None, help='database of contaminant '
                             'sequences in Fasta/Fastq format')
    contam_args.add_argument('--contam-memory', metavar='MEM', default='1e6',
                             type=khmer_args.memory_setting,
                             help='memory to allocate for storing contaminant '
                             'sequences; default is 1M')
    contam_args.add_argument('--contam-max-fpr', type=float, metavar='FPR',
                             default=0.001, help='terminate if the expected '
                             'false positive rate is higher than the specified'
                             ' FPR; default is 0.001')

    filter_args = subparser.add_argument_group(
        'Filtering k-mers',
        'Memory constraints often require running `kevlar novel` with false '
        'positive rates (FPRs) in the 0.1 - 0.2 range, resulting in some '
        'k-mers reported with highly inflated abundances. This script handles '
        'a much smaller amount of data, and in limited memory can achieve a '
        'much lower FPR, compute exact k-mer abundances, and discard '
        '"interesting" k-mers whose abundances were incorrectly reported '
        'previously.'
    )
    filter_args.add_argument('--abund-memory', metavar='MEM', default='1e6',
                             type=khmer_args.memory_setting,
                             help='memory to allocate for re-calculating '
                             'abundance of interesting k-mers; default is 1M')
    filter_args.add_argument('--abund-max-fpr', type=float, metavar='FPR',
                             default=0.001, help='terminate if the expected '
                             'false positive rate is higher than the specified'
                             ' FPR; default is 0.001')
    filter_args.add_argument('--min-abund', type=int, default=5, metavar='Y',
                             help='minimum abundance required to call a '
                             'k-mer novel; should be the same value used for '
                             '--case_min in `kevlar novel`; default is 5')
    filter_args.add_argument('--skip2', default=False, action='store_true',
                             help='skip the second pass over the reads that '
                             'recalculates abundance after reference and '
                             'contaminant k-mers are discarded')
    filter_args.add_argument('--ignore', metavar='KM', nargs='+',
                             help='ignore the specified k-mer(s)')

    misc_args = subparser.add_argument_group(
        'Miscellaneous settings'
    )
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    misc_args.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    misc_args.add_argument('--aug-out', metavar='FILE',
                           help='optional augmented Fastq output')
    misc_args.add_argument('--cc-prefix', metavar='PREFIX',
                           help='group reads by novel k-mers, and use the '
                           'specified prefix to write each group to its own '
                           'file')

    subparser.add_argument('augfastq', nargs='+', help='one or more files in '
                           '"augmented" Fastq format (a la `kevlar novel` '
                           'output)')


def load_refr(refrfile, ksize, memory, maxfpr=0.001, logfile=sys.stderr):
    """Load reference genome or contaminant database from a file."""
    buckets = memory * khmer._buckets_per_byte['nodegraph'] / 4
    refr = khmer.Nodetable(ksize, buckets, 4)
    nr, nk = refr.consume_seqfile(refrfile)
    fpr = kevlar.calc_fpr(refr)
    message = '    {:d} sequences and {:d} k-mers consumed'.format(nr, nk)
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    print(message, file=logfile)
    if fpr > maxfpr:
        print('[kevlar::filter] FPR too high, bailing out', file=logfile)
        sys.exit(1)
    return refr


def load_input(filelist, ksize, memory, maxfpr=0.001, logfile=sys.stderr):
    """
    Load input data.

    The input data is loaded into two data structures. First, the read
    sequences are loaded into a countgraph to recompute k-mer abundances with
    (effectively) exact precision. Second, the reads and their corresponding
    "interesting" k-mers are loaded into an AnnotatedReadSet to de-duplicate
    reads and group k-mers by read.
    """
    countgraph = khmer.Countgraph(ksize, memory / 4, 4)
    read_inst_consumed = 0
    int_kmer_instances = 0
    int_kmers_parsed = set()
    readset = kevlar.seqio.AnnotatedReadSet()
    for filename in filelist:
        print('    -', filename, file=logfile)
        with kevlar.open(filename, 'r') as infile:
            for record in kevlar.parse_augmented_fastq(infile):
                if record.name not in readset._reads:
                    countgraph.consume(record.sequence)
                readset.add(record)
                read_inst_consumed += 1
                for kmer in record.ikmers:
                    int_kmer_instances += 1
                    minkmer = kevlar.revcommin(kmer.sequence)
                    int_kmers_parsed.add(minkmer)
    n_kmers_distinct = len(int_kmers_parsed)

    fpr = kevlar.calc_fpr(countgraph)
    message = '    {:d} instances'.format(read_inst_consumed)
    message += ' of {:d} reads consumed'.format(len(readset))
    message += ', annotated with {:d} instances '.format(int_kmer_instances)
    message += 'of {:d} distinct "interesting" k-mers'.format(n_kmers_distinct)
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    print(message, file=logfile)
    if fpr > maxfpr:
        print('[kevlar::filter] FPR too high, bailing out', file=logfile)
        sys.exit(1)
    return readset, countgraph


def validate_and_print(readset, countgraph, refr=None, contam=None, minabund=5,
                       skip2=False, outfile=sys.stdout, augout=None,
                       logfile=sys.stderr):
    readset.validate(countgraph, refr=refr, contam=contam, minabund=minabund)
    if not skip2:
        ksize, tablesizes = countgraph.ksize(), countgraph.hashsizes()
        countgraph = khmer._Countgraph(ksize, tablesizes)
        readset.recalc_abund(countgraph, minabund)

    n = 0  # Get an unbound var error later (printing report) without this?!?!
    for n, record in enumerate(readset):
        khmer.utils.write_record(record, outfile)
        if augout:
            kevlar.print_augmented_fastq(record, augout)

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
    message += '{:d} contaminant reads discarded'.format(readset.contam)
    message += '\n        '
    message += '{:d} reads written to output'.format(n + 1)
    print(message, file=logfile)


def main(args):
    if args.cc_prefix:  # pragma: no cover
        try:
            import networkx
        except ImportError:
            print('[kevlar::filter] FATAL ERROR: cannot group reads by novel '
                  'k-mers (--cc-prefix flag) unless the "networkx" module is '
                  'installed', file=sys.stderr)
            sys.exit(1)

    timer = kevlar.Timer()
    timer.start()

    refr = None
    if args.refr:
        timer.start('loadrefr')
        print('[kevlar::filter] Loading reference genome from',
              args.refr, file=args.logfile)
        refr = load_refr(args.refr, args.ksize, args.refr_memory,
                         args.refr_max_fpr, args.logfile)
        elapsed = timer.stop('loadrefr')
        print('[kevlar::filter]',
              'Reference genome loaded in {:.2f} sec'.format(elapsed),
              file=args.logfile)

    contam = None
    if args.contam:
        timer.start('loadcontam')
        print('[kevlar::filter] Loading contaminants from', args.contam,
              file=args.logfile)
        contam = load_refr(args.contam, args.ksize, args.contam_memory,
                           args.contam_max_fpr, args.logfile)
        elapsed = timer.stop('loadcontam')
        print('[kevlar::filter]',
              'Contaminant database loaded in {:.2f} sec'.format(elapsed),
              file=args.logfile)

    timer.start('recalc')
    print('[kevlar::filter] Loading input; recalculate k-mer abundances,',
          'de-duplicate reads and merge k-mers',
          file=args.logfile)
    readset, countgraph = load_input(args.augfastq, args.ksize,
                                     args.abund_memory, args.abund_max_fpr,
                                     args.logfile)
    elapsed = timer.stop('recalc')
    print('[kevlar::filter] Input loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    timer.start('validate')
    print('[kevlar::filter] Validate k-mers and print reads',
          file=args.logfile)
    outstream = kevlar.open(args.out, 'w')
    augstream = kevlar.open(args.aug_out, 'w') if args.aug_out else None
    validate_and_print(readset, countgraph, refr, contam, args.min_abund,
                       args.skip2, outstream, augstream, args.logfile)
    elapsed = timer.stop('validate')
    print('[kevlar::filter] k-mers validated and reads printed',
          'in {:.2f} sec'.format(elapsed), file=args.logfile)

    if args.cc_prefix:
        timer.start('graph')
        print('[kevlar::filter] Group reads by novel k-mers',
              file=args.logfile)
        readset.group_reads_by_novel_kmers(args.cc_prefix,
                                           logstream=args.logfile)
        elapsed = timer.stop('graph')
        print('[kevlar::filter] reads grouped by novel k-mers',
              'in {:.2f} sec'.format(elapsed), file=args.logfile)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::filter]', message, file=args.logfile)
