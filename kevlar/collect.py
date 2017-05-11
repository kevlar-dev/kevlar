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
    subparser = subparsers.add_parser('collect')
    subparser.add_argument('-d', '--debug', action='store_true',
                           help='print debugging output')
    subparser.add_argument('-M', '--memory', default='1e6', metavar='MEM',
                           type=khmer_args.memory_setting,
                           help='memory to allocate for recalculating '
                           'abundances of novel k-mers; default is 1M')
    subparser.add_argument('--minabund', type=int, default=5, metavar='Y',
                           help='minimum case abundance required to call a '
                           'k-mer novel; used to filter out k-mers with '
                           'inflated abundances; should equal the value of '
                           '--case_min used in "kevlar novel"; default is 5')
    subparser.add_argument('--max-fpr', type=float, metavar='FPR',
                           default=0.001, help='terminate if the expected '
                           'false positive rate is higher than the specified '
                           'FPR; default is 0.001')
    subparser.add_argument('--refr', metavar='FILE', type=str, default=None,
                           help='reference genome in Fasta/Fastq format; any '
                           'k-mers designated as "interesting" by "kevlar '
                           'novel" are ignored if they are present in the '
                           'reference genome')
    subparser.add_argument('--refr-memory', metavar='MEM', default='1e6',
                           type=khmer_args.memory_setting,
                           help='memory to allocate for storing the reference '
                           'genome; default is 1M')
    subparser.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    subparser.add_argument('--ignore', metavar='KMER', nargs='+',
                           help='ignore the specified k-mer(s)')
    subparser.add_argument('-o', '--out', type=argparse.FileType('w'),
                           metavar='OUT',
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('--collapse', action='store_true', help='collapse '
                           'linear paths contained in other linear paths')
    subparser.add_argument('novel_output', nargs='+', help='one or more output'
                           ' files from the "kevlar novel" command')


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
        with kevlar.open(filename, 'r') as infile:
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
        print('[kevlar::collect] FPR too high, bailing out', fpr, file=logfile)
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
        with kevlar.open(filename, 'r') as infile:
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
    if args.refr:
        refr = load_refr(args.refr, args.ksize, args.refr_memory)
    countgraph = recalc_abund(args.novel_output, args.ksize, args.memory,
                              args.max_fpr, args.logfile)
    variants = load_novel_kmers(args.novel_output, countgraph, refr,
                                args.minabund, args.logfile)
    assemble_contigs(countgraph, variants, args.ignore, args.collapse,
                     args.debug, args.logfile)
    variants.write(outstream=args.out)
