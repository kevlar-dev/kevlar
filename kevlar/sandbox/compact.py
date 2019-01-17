#!/usr/bin/env python
# -*- coding: utf-8 -*-
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


cli = argparse.ArgumentParser()
cli.add_argument('-d', '--delta', metavar='Δ', type=int, default=10, help='support approximate matches by extending each query by Δ bp each direction')
cli.add_argument('reference', help='location of reference calls in BED format')
cli.add_argument('vcf', help='variant predictions in VCF format')
args = cli.parse_args()

bedstream = kevlar.open(args.reference, 'r')
index = kevlar.evaluate.populate_index_from_bed(bedstream)

vcfstream = kevlar.open(args.vcf, 'r')
reader = kevlar.vcf.VCFReader(vcfstream)
calls = list(reader)

writer = kevlar.vcf.VCFWriter(sys.stdout, source='kevlar::sandbox::compact.py')
writer.register_samples_from_reader(reader)
writer.write_header()
for varcall in kevlar.evaluate.compact(calls, index, delta=args.delta):
    writer.write(varcall)
