#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019 Battelle National Biodefense Institute
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
import kevlar
import khmer
import sys


allocators = {
    'nodegraph': khmer.Nodegraph,
    'countgraph': khmer.Countgraph,
    'smallcountgraph': khmer.SmallCountgraph,
    'nodetable': khmer.Nodetable,
    'counttable': khmer.Counttable,
    'smallcounttable': khmer.SmallCounttable,
}


cli = argparse.ArgumentParser()
cli.add_argument('--sketch-type', metavar='T', choices=allocators.keys(), default='counttable', help='Sketch type to use for output')
cli.add_argument('--num-tables', type=int, default=4, metavar='N')
cli.add_argument('--table-size', type=int, default=1000, metavar='X')
cli.add_argument('sketch', help='original sketch')
cli.add_argument('subsketch', help='new sketch to create')
cli.add_argument('sequence', nargs='+', help='sequences to sample from sketch')
args = cli.parse_args()

sketch = kevlar.sketch.load(args.sketch)
allocfunc = allocators[args.sketch_type]
subsketch = allocfunc(sketch.ksize(), args.table_size, args.num_tables)

kmers = set()
for seq in args.sequence:
    for kmer in sketch.get_kmers(seq):
        minkmer = kevlar.revcommin(kmer)
        kmers.add(minkmer)
for kmer in kmers:
    count = sketch.get(kmer)
    for _ in range(count):
        subsketch.add(kmer)

subsketch.save(args.subsketch)
fpr = khmer.calc_expected_collisions(subsketch, max_false_pos=100.0)
print('Estimated FPR: {:.4f}'.format(fpr))
