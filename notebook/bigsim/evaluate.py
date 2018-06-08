#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
from collections import defaultdict
import sys

import intervaltree
from intervaltree import IntervalTree
import pandas
from evalutils import IntervalForest, populate_index_from_simulation, compact
from evalutils import assess_variants, assess_variants_mvf
import kevlar
from kevlar.vcf import VCFReader


parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tolerance', type=int, metavar='T', default=10,
                    help='extend real variants by T nucleotides when querying '
                    'for overlap with variant calls; default is 10')
parser.add_argument('--mvf', action='store_true', help='input is in .mvf '
                    'format; default is .vcf')
parser.add_argument('--correct', help='print correct variants to file')
parser.add_argument('--missing', help='print missing variants to file')
parser.add_argument('--false', help='print false variants to file')
parser.add_argument('--collisions', help='print calls that match the same '
                    'variant')
parser.add_argument('simvar', help='simulated variants (in custom 3-4 column '
                    'tabular format)')
parser.add_argument('varcalls', help='VCF file of variant calls')
args = parser.parse_args()


index = populate_index_from_simulation(args.simvar, 'chr17')
if args.mvf:
    table = pandas.read_table(args.varcalls, sep='\t')
    variants = table.sort_values('CHILD_DP', ascending=False)
    assess_func = assess_variants_mvf
else:
    reader = VCFReader(kevlar.open(args.varcalls, 'r'))
    variants = compact(reader, index, delta=args.tolerance)
    assess_func = assess_variants
correct, false, missing, mapping = assess_func(
    variants, index, delta=args.tolerance
)


numcollisions = 0
for variant, calllist in mapping.items():
    if len(calllist) > 1:
        numcollisions += 1
if numcollisions > 0:
    print('WARNING:', numcollisions, 'variants matched by multiple calls',
          file=sys.stderr)
    if args.collisions:
        with open(args.collisions, 'w') as outstream:
            for variant, calllist in mapping.items():
                if len(calllist) > 1:
                    print('\n#VARIANT:', variant, file=outstream)
                    for varcall in calllist:
                        if args.mvf:
                            print('    -', varcall, file=outstream)
                        else:
                            print('    -', varcall.vcf, file=outstream)

if args.missing:
    with open(args.missing, 'w') as outstream:
        for variant in missing:
            print(variant.begin, *variant.data.split('<-'), sep='\t',
                  file=outstream)

if args.correct:
    outstream = kevlar.open(args.correct, 'w')
    if args.mvf:
        for varcall in correct:
            print(varcall, file=outstream)
    else:
        writer = kevlar.vcf.VCFWriter(outstream)
        for varcall in correct:
            writer.write(varcall)

if args.false:
    outstream = kevlar.open(args.false, 'w')
    if args.mvf:
        for varcall in false:
            print(varcall, file=outstream)
    else:
        writer = kevlar.vcf.VCFWriter(outstream)
        for varcall in false:
            writer.write(varcall)

print('Correct:', len(mapping))
print('False:', len(false))
print('Missing:', len(missing))
