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


class ReadSet(object):

    def __init__(self):
        self._reads = dict()
        self._kmers = dict()

        self._masked_count = 0
        self._masked_distinct = 0
        self._lowabund_count = 0
        self._lowabund_distinct = 0
        self._valid_count = 0
        self._valid_distinct = 0
        self._novalidkmers_count = 0

    def __len__(self):
        assert len(self._reads) == len(self._kmers)
        return len(self._reads)

    @property
    def interesting(self):
        instances = sum([
            self._masked_count,
            self._lowabund_count,
            self._valid_count
        ])
        distinct = sum([
            self._masked_distinct,
            self._lowabund_distinct,
            self._valid_distinct
        ])
        return instances, distinct

    def __iter__(self):
        for readid in self._reads:
            if len(self._kmers[readid]) == 0:
                self._novalidkmers_count += 1
            else:
                yield self._reads[readid], self._kmers[readid]

    def add(self, record, kmers):
        if record.name in self._reads:
            assert record.sequence == self._reads[record.name].sequence
            self._kmers[record.name].update(kmers)
        else:
            self._reads[record.name] = record
            self._kmers[record.name] = kmers

    def validate(self, counts, refr=None, minabund=5):
        masked = set()
        lowabund = set()
        valid = set()
        for readid in self._kmers:
            validated_kmers = dict()
            for offset in self._kmers[readid]:
                kmer, abundances = self._kmers[readid][offset]
                if refr and refr.get(kmer) > 0:
                    self._masked_count += 1
                    if kmer not in masked:
                        self._masked_distinct += 1
                        masked.add(kmer)
                elif counts.get(kmer) < minabund:
                    self._lowabund_count += 1
                    if kmer not in lowabund:
                        self._lowabund_distinct += 1
                        lowabund.add(kmer)
                else:
                    validated_kmers[offset] = (kmer, abundances)
                    self._valid_count += 1
                    if kmer not in valid:
                        self._valid_distinct += 1
                        valid.add(kmer)
            self._kmers[readid] = validated_kmers


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

    filter_args = subparser.add_argument_group(
        'Filtering k-mers',
        'Memory constraints often require running `kevlar find` with false '
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
                             '--case_min in `kevlar find`; default is 5')
    filter_args.add_argument('--ignore', metavar='KM', nargs='+',
                             help='ignore the specified k-mer(s)')

    misc_args = subparser.add_argument_group(
        'Miscellaneous settings'
    )
    misc_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
    misc_args.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    misc_args.add_argument('-o', '--out', type=argparse.FileType('w'),
                           metavar='FILE',
                           help='output file; default is terminal (stdout)')
    misc_args.add_argument('--aug-out', type=argparse.FileType('w'),
                           metavar='FILE',
                           help='optional augmented Fastq output')

    subparser.add_argument('augfastq', nargs='+', help='one or more files in '
                           '"augmented" Fastq format (a la `kevlar find` '
                           'output)')


def load_refr(refrfile, ksize, memory, maxfpr=0.001, logfile=sys.stderr):
    buckets = memory * khmer._buckets_per_byte['nodegraph'] / 4
    refr = khmer.Nodetable(ksize, int(buckets), 4)
    nr, nk = refr.consume_seqfile(refrfile)
    fpr = kevlar.calc_fpr(refr)
    message = '    {:d} sequences and {:d} k-mers consumed'.format(nr, nk)
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    print(message, file=logfile)
    if fpr > maxfpr:
        print('[kevlar::filter] FPR too high, bailing out', file=logfile)
        sys.exit(1)
    return refr


def recalc_abund(filelist, ksize, memory, maxfpr=0.001, logfile=sys.stderr):
    countgraph = khmer.Countgraph(ksize, memory / 4, 4)
    read_inst_consumed = 0
    int_kmer_parsed = 0
    readset = ReadSet()
    for filename in filelist:
        print('    -', filename, file=logfile)
        with open(filename, 'r') as infile:
            for record, kmers in kevlar.parse_augmented_fastq(infile):
                if record.name not in readset._reads:
                    countgraph.consume(record.sequence)
                readset.add(record, kmers)
                read_inst_consumed += 1

    fpr = kevlar.calc_fpr(countgraph)
    message = '    {:d} instances'.format(read_inst_consumed)
    message += ' of {:d} reads consumed'.format(len(readset))
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    print(message, file=logfile)
    if fpr > maxfpr:
        print('[kevlar::filter] FPR too high, bailing out', file=logfile)
        sys.exit(1)
    return readset, countgraph


def validate_and_print(readset, countgraph, refr=None, minabund=5,
                       outfile=sys.stdout, augout=None, logfile=sys.stderr):
    readset.validate(countgraph, refr, minabund)
    for record, kmers in readset:
        khmer.utils.write_record(record, outfile)
        if augout:
            khmer.utils.write_record(record, augout)
            for offset in sorted(kmers):
                kmer, abundances = kmers[offset]
                abundstr = ' '.join([str(a) for a in abundances])
                print(' ' * offset, kmer, ' ' * 10, abundstr, '#', sep='',
                      file=augout)

    int_instances, int_distinct = readset.interesting
    message = '    processed {:d} instances'.format(int_instances)
    message += ' of {:d} distinct "interesting" k-mers'.format(int_distinct)
    message += ' in {:d} reads'.format(len(readset))
    message += '\n        '
    message += '{:d} instances'.format(readset._masked_count)
    message += ' of {:d} distinct k-mers'.format(readset._masked_distinct)
    message += ' masked by the reference genome'
    message += '\n        '
    message += '{:d} instances'.format(readset._lowabund_count)
    message += ' of {:d} distinct k-mers'.format(readset._lowabund_distinct)
    message += ' discarded due to low abundance'
    message += '\n        '
    message += '{:d} instances'.format(readset._valid_count)
    message += ' of {:d} distinct k-mers'.format(readset._valid_distinct)
    message += ' validated as novel'
    message += '\n        '
    message += '{:d} reads'.format(readset._novalidkmers_count)
    message += ' with no surviving valid k-mers ignored'
    print(message, file=logfile)


def main(args):
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

    timer.start('recalc')
    print('[kevlar::filter] Loading input; recalculate k-mer abundances',
          'de-duplicate reads and merge k-mers',
          file=args.logfile)
    readset, countgraph = recalc_abund(args.augfastq, args.ksize,
                                       args.abund_memory, args.abund_max_fpr,
                                       args.logfile)
    elapsed = timer.stop('recalc')
    print('[kevlar::filter] Input loaded in {:.2f} sec'.format(elapsed))

    timer.start('validate')
    print('[kevlar::filter] Validate k-mers and print reads',
          file=args.logfile)
    validate_and_print(readset, countgraph, refr, args.min_abund, args.out,
                       args.aug_out, args.logfile)
    elapsed = timer.stop('validate')
    print('[kevlar::filter] k-mers validated and reads printed',
          'in {:.2f} sec'.format(elapsed), file=args.logfile)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::filter]', message, file=args.logfile)
