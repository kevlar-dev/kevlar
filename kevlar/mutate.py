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

char_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
index_to_char = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
mutation_functions = {
    'snv': mutate_snv,
    'ins': mutate_insertion,
    'del': mutate_deletion,
    'inv': mutate_inversion,
}


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


def mutate_snv(sequence, mutation):
    refrbase = sequence[mutation.pos]
    nuclindex = char_to_index[refrbase]
    newindex = nuclindex + int(mutation.data)
    while newindex > 3:
        newindex -= 4
    newbase = index_to_char[newindex]
    prefix, suffix = sequence[:mutation.pos], sequence[mutation.pos+1:]
    return prefix + newbase + suffix


def mutate_insertion(sequence, mutation):
    prefix, suffix = sequence[:mutation.pos], sequence[mutation.pos:]
    return prefix + mutation.data + suffix


def mutate_deletion(sequence, mutation):
    del_length = int(mutation.data)
    prefix = sequence[:mutation.pos]
    suffix = sequence[mutation.pos+del_length:]
    return prefix + suffix


def mutate_inversion(sequence, mutation):
    inv_length = int(mutation.data)
    prefix = sequence[:mutation.pos]
    suffix = sequence[mutation.pos+del_length:]
    invseq = sequence[mutation.pos:mutation.pos+inv_length:-1]
    return prefix + invseq + suffix


def mutate_sequence(sequence, mutlist):
    for mutation in mutlist:
        mutfunc = mutation_functions[mutation.type]
        sequence = mutfunc(sequence, mutation)
    return sequence


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
