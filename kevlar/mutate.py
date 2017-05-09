#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict, NamedTuple
import argparse
import sys
import khmer
from khmer.utils import write_record
import kevlar


Mutation = NamedTuple('Mutation', 'seq pos type data')


# SNV: chr1 5087 snv 2
# insertion: scaffold14 19983 ins AGCTACCCCAGTGAGTCGGTAATGTGATC
# deletion: contig8 8837 del 5
# inversion: X 2884322 inv 2766


def subparser(subparsers):
    subparser = subparsers.add_parser('mutate')
    subparser.add_argument('-o', '--out', metavar='FILE', default=sys.stdout,
                           type=argparse.FileType('w'),
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('mutations', type=argparse.FileType('r'),
                           help='mutations file')
    subparser.add_argument('genome', metavar='FILE', help='genome to mutate')


def load_mutations(instream, logstream):
    mutations = defaultdict(list)
    count = 0
    for line in instream:
        if line.startswith('#') or line.strip() == '':
            continue
        try:
            sequence, offset, vartype, data = line.strip().split()
        except ValueError:
            raise ValueError('error parsing mutation: ' + line)
        if vartype not in ['snv', 'ins', 'del', 'inv']:
            raise ValueError('invalid variant type "{:s}"'.format(vartype))
        mut = Mutation(seq=sequence, pos=int(offset), type=vartype, data=data)
        mutations[sequence].append(mut)
        count += 1
    message = '    loaded {:d} mutations'.format(count)
    message += ' on {:d} sequences'.format(len(mutations), file=logstream)
    return mutations


def mutate_sequence(sequence, mutlist):
    for mutation in mutlist:
        if mutation.type == 'snv':
            pass


def mutate_genome(infile, mutations):
    parser = khmer.ReadParser(infile)
    for record in parser:
        if record.name in mutations:
            mutlist = sorted(mutations[record.name], key=lambda m: m.pos,
                             reverse=True)
            record.sequence = mutate_sequence(record.sequence, mutlist)
        yield record


def main(args):
    print('[kevlar::mutate] loading mutations', file=args.logfile)
    mutations = load_mutations(args.mutations, args.logstream)

    print('[kevlar::mutate] mutating genome', file=args.logfile)
    for record in mutate_genome(args.genome, mutations):
        write_record(record, args.out)
