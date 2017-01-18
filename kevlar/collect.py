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
    subparser.add_argument('-M', '--memory', default='1e6',
                           type=khmer_args.memory_setting, metavar='MEM',
                           help='total memory to allocate for the node '
                           'graph; default is 1M')
    subparser.add_argument('--max-fpr', type=float, default=0.2, metavar='FPR',
                           help='terminate if the expected false positive rate'
                           ' is higher than the specified FPR; default is 0.2')
    subparser.add_argument('-k', '--ksize', type=int, default=31, metavar='K',
                           help='k-mer size; default is 31')
    subparser.add_argument('--ignore', metavar='KMER', nargs='+',
                           help='ignore the specified k-mer(s)')
    subparser.add_argument('-o', '--out', type=argparse.FileType('w'),
                           metavar='OUT',
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('--collapse', action='store_true', help='collapse '
                           'linear paths contained in other linear paths')
    subparser.add_argument('find_output', nargs='+',
                           type=argparse.FileType('r'), help='one or more '
                           'output files from the "kevlar find" command')


def load_input(infile, nodegraph, variants):
    readid = None
    reads_consumed = 0
    kmers_consumed = 0
    for line in infile:
        if line.startswith('@'):
            # Load the 4 lines of the Fastq record
            readid = line[1:].strip()
            seq = next(infile).strip()
            _ = next(infile)
            _ = next(infile)

            nkmers = nodegraph.consume(seq)
            reads_consumed += 1
            kmers_consumed += nkmers
        elif line.endswith('#\n'):
            novel_kmer = re.search('^ *(\S+)', line).group(1)
            variants.add_kmer(novel_kmer, readid)

    return reads_consumed, kmers_consumed


def load_all_inputs(filelist, nodegraph, variants, maxfpr=0.2,
                    logfile=sys.stderr):
    kmers_consumed = 0
    reads_consumed = 0
    for infile in filelist:
        infilename = 'input file'
        if hasattr(infile, "name"):
            infilename = infile.name
        print('[kevlar::collect] Loading', infilename, file=logfile)
        nr, nk = load_input(infile, nodegraph, variants)
        reads_consumed += nr
        kmers_consumed += nk
    fpr = kevlar.calc_fpr(nodegraph)
    message = '    {:d} reads'.format(reads_consumed)
    message += ' and {:d} k-mers consumed'.format(kmers_consumed)
    message += '; {:d} instances'.format(variants.nkmerinst)
    message += ' of {:d} unique novel k-mers'.format(variants.nkmers)
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    print(message, file=logfile)
    if fpr > maxfpr:
        sys.exit(1)


def assemble_contigs(nodegraph, variants, kmers_to_ignore=None,
                     collapse=True, debug=False, logfile=sys.stderr):
    print('[kevlar::collect] Retrieving linear paths', file=logfile)
    for kmer in variants.kmers:
        if kmers_to_ignore and kmer in kmers_to_ignore:
            continue
        if debug:
            print('[kevlar::collect]     DEBUG kmer:', kmer, file=logfile)
        contig = nodegraph.assemble_linear_path(kmer)
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
    nodegraph = khmer.Nodegraph(args.ksize, args.memory / 4, 4)
    variants = kevlar.VariantSet()

    load_all_inputs(args.find_output, nodegraph, variants, args.max_fpr,
                    args.logfile)
    assemble_contigs(nodegraph, variants, args.ignore, args.collapse,
                     args.debug, args.logfile)
    variants.write(outstream=args.out)
