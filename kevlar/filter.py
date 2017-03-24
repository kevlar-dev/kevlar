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
        'A k-mer labeled as "interesting" by `kevlar find` should be discarded'
        ' if it is present in the supplied reference genome.'
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
        'positive rates (FPRs) in the 0.1 - 0.2 range. This will result in '
        'some k-mers with highly inflated abundances. With much less input '
        'data, this script can achieve a much lower FPR --> exact k-mer '
        'abundances.'
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
                           metavar='OUT',
                           help='output file; default is terminal (stdout)')

    subparser.add_argument('find_output', nargs='+', help='one or more output '
                           'files from the `kevlar find` command')


def load_refr(refrfile, ksize, memory, logfile=sys.stderr):
    print('[kevlar::collect] Loading reference from ' + refrfile, file=logfile)
    buckets = memory * khmer._buckets_per_byte['nodegraph'] / 4
    refr = khmer.Nodetable(ksize, int(buckets), 4)
    nr, nk = refr.consume_seqfile(refrfile)
    message = '    {:d} reads and {:d} k-mers consumed'.format(nr, nk)
    print(message, file=logfile)
    return refr


def recalc_abund(filelist, ksize, memory, maxfpr=0.1, logfile=sys.stderr):
    countgraph = khmer.Countgraph(ksize, memory / 4, 4)
    kmers_consumed = 0
    reads_consumed = 0
    readids = set()
    print('[kevlar::collect]',
          'Loading input, pass 1; recalculating k-mer abundances',
          file=logfile)
    for filename in filelist:
        print('   ', filename, file=logfile)
        with open(filename, 'r') as infile:
            for line in infile:
                if line.startswith('@'):
                    # Load the 4 lines of the Fastq record
                    readid = line[1:].strip()
                    seq = next(infile).strip()
                    _ = next(infile)
                    _ = next(infile)

                    if readid not in readids:
                        readids.add(readid)
                        kmers_consumed += countgraph.consume(seq)
                        reads_consumed += 1

                elif line.endswith('#\n'):
                    # Ignore the novel k-mers in the first pass
                    continue

    fpr = kevlar.calc_fpr(countgraph)
    message = '    {:d} reads'.format(reads_consumed)
    message += ' and {:d} k-mers consumed'.format(kmers_consumed)
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    print(message, file=logfile)
    if fpr > maxfpr:
        sys.exit(1)
    return countgraph


def load_novel_kmers(filelist, countgraph, refr=None, minabund=5,
                     logfile=sys.stderr):
    variants = kevlar.VariantSet()
    print('[kevlar::collect] Loading input, pass 2; store novel k-mers',
          file=logfile)
    masked = 0
    discarded = 0
    for filename in filelist:
        print('   ', filename, file=logfile)
        with open(filename, 'r') as infile:
            for line in infile:
                if line.startswith('@'):
                    # Load the 4 lines of the Fastq record
                    readid = line[1:].strip()
                    seq = next(infile).strip()
                    _ = next(infile)
                    _ = next(infile)
                elif line.endswith('#\n'):
                    # Second pass, double check the abundance, and then add
                    novel_kmer = str(re.search('^ *(\S+)', line).group(1))
                    if refr and refr.get(novel_kmer) > 0:
                        masked += 1
                    elif countgraph.get(novel_kmer) < minabund:
                        discarded += 1
                    else:
                        variants.add_kmer(novel_kmer, readid)

    message = '    found {:d} instances'.format(variants.nkmerinst)
    message += ' of {:d} unique novel k-mers'.format(variants.nkmers)
    message += '; {:d} k-mers masked by the reference genome'.format(masked)
    message += '; {:d} k-mers with inflated counts discarded'.format(discarded)
    print(message, file=logfile)
    return variants


def assemble_contigs(countgraph, variants, kmers_to_ignore=None,
                     collapse=True, debug=False, logfile=sys.stderr):
    asm = khmer.JunctionCountAssembler(countgraph)
    print('[kevlar::collect] Retrieving linear paths', file=logfile)
    for kmer in variants.kmers:
        if kmers_to_ignore and kmer in kmers_to_ignore:
            continue
        if debug:
            print('[kevlar::collect]     DEBUG kmer:', kmer, file=logfile)
        contigs = asm.assemble(kmer)
        for i, contig in enumerate(contigs):
            if contig == '':
                print('    WARNING: no linear path found for k-mer', kmer,
                      file=logfile)
                continue
            variants.add_contig(contig, kmer)
    print('    {:d} linear paths'.format(variants.ncontigs), file=logfile)

    if collapse:
        print('[kevlar::collect] Collapsing linear paths', file=logfile)
        variants.collapse()
        print('    {:d} collapsed linear paths'.format(variants.ncontigs),
              file=logfile)


def main(args):
    refr = None
    if args.mask:
        refr = load_refr(args.refr, args.ksize, args.refr_memory)
    countgraph = recalc_abund(args.find_output, args.ksize, args.memory,
                              args.max_fpr, args.logfile)
    variants = load_novel_kmers(args.find_output, countgraph, refr,
                                args.minabund, args.logfile)
    assemble_contigs(countgraph, variants, args.ignore, args.collapse,
                     args.debug, args.logfile)
    variants.write(outstream=args.out)
