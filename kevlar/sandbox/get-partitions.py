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
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--out', metavar='FILE', help='output filename')
parser.add_argument('-p', '--out-pattern', metavar='REGEX', help='out file name pattern with a {} placeholder for partition ID')
parser.add_argument('augfastx')
parser.add_argument('partition', nargs='+')
args = parser.parse_args()

if args.out and args.out_pattern:
    raise Exception('cannot give outfile and outpattern together')
elif args.out and not args.out_pattern:
    args.out = kevlar.open(args.out, 'w')
elif not args.out and not args.out_pattern:
    args.out = sys.stdout

partids = set(args.partition)
fh = kevlar.open(args.augfastx, 'r')
reader = kevlar.parse_augmented_fastx(fh)
preader = kevlar.parse_partitioned_reads(reader)
for partid, partition in preader:
    if partid not in partids:
        continue
    if args.out_pattern:
        pattern = str(args.out_pattern)
        outfile = pattern.format(partid)
        with kevlar.open(outfile, 'w') as out:
            for read in partition:
                kevlar.print_augmented_fastx(read, out)
    else:
        for read in partition:
            kevlar.print_augmented_fastx(read, args.out)
