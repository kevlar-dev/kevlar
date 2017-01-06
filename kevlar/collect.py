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
import sys

import khmer
from khmer import khmer_args
import kevlar


def subparser(subparsers):
    subparser = subparsers.add_parser('collect')
    subparser.add_argument('-M', '--memory', default='1e6',
                           type=khmer_args.memory_setting, metavar='MEM',
                           help='total memory to allocate for the node '
                           'graph; default is 1M')
    subparser.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    subparser.add_argument('--out', type=argparse.FileType('w'),
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('--collapse', action='store_true', help='collapse '
                           'linear paths contained in other linear paths')
    subparser.add_argument('find_output', nargs='+', help='one or more output '
                           'files from the "kevlar find" command')


def load_input(infile, nodegraph, unique_kmers):
    readid = None
    reads_consumed = 0
    kmers_consumed = 0
    kmer_instances = 0
    for line in infile:
        if line.startswith('@'):
            # Load the 4 lines of the Fastq record
            readid = line[1:].strip()
            seq = next(infile).strip()
            _ = next(infile)
            _ = next(infile)

            _, nkmers = nodegraph.consume(seq)
            reads_consumed += 1
            kmers_consumed += nkmers
        elif line.endswith('#\n'):
            novel_kmer = re.search('^ *(\S+)', line)
            kmer_instances += 1
            minkmer = kevlar.revcommin(novel_kmer)
            unique_kmers.add(minkmer)
            kmer_reads[minkmer].append(readid)

    return kmers_consumed, reads_consumed, kmer_instances, unique_kmers


def main(args):
    nodegraph = khmer.Nodegraph(args.ksize, args.memory / 4, 4)

    kmers_consumed = 0
    reads_consumed = 0
    novel_kmer_instances = 0
    unique_novel_kmers = set()
    kmer_reads = defaultdict(list)
    for infile in args.find_output:
        print('[kevlar::collect] Loading', infile, file=args.logfile)
        with open(infile, 'r') as fh:
            data = load_input(fh, nodegraph, unique_novel_kmers, kmer_reads)
            kmers_consumed += data[0]
            reads_consumed += data[1]
            novel_kmer_instances += data[2]
    message = '    {:d} reads'.format(reads_consumed)
    message += ' and {:d} k-mers consumed'.format(kmers_consumed)
    message += '; {:d} instances'.format(novel_kmer_instances)
    message += ' of {:d} unique novel k-mers'.format(len(unique_novel_kmers))
    print(message, file=args.logfile)

    print('[kevlar::collect] Retrieving linear paths', file=args.logfile)
    variants = kevlar.VariantSet()
    for kmer in unique_novel_kmers:
        linear_path = nodegraph.assemble_linear_path(kmer)
        if linear_path == '':
            print('    WARNING: no linear path found for k-mer', kmer,
                  file=args.logfile)
            continue
        for readid in kmer_reads[kmer]:
            variants.add_kmer(kmer, readid, linear_path)
    print('    {:d} linear paths'.format(variants.npaths), file=args.logfile)

    if args.collapse:
        print('[kevlar::collect] Collapsing linear paths', file=args.logfile)
        variants.collapse()
        print('    {:d} collapsed linear paths'.format(variants.npaths),
              file=args.logfile)

    variants.path_table(outstream=args.out)
