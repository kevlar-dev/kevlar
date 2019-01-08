#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019 Battelle National Biodefense Institute
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar


def load_variant_mask(bedstream):
    index = kevlar.IntervalForest()
    for chrom, start, end, data in kevlar.parse_bed(bedstream):
        index.insert(chrom, start, end)
    return index


def varfilter(calls, varmask):
    for varcall in calls:
        if varmask.query(varcall.seqid, varcall.position) != set():
            varcall.filter(kevlar.vcf.VariantFilter.UserFilter)
        yield varcall


def main(args):
    reader = kevlar.vcf.vcfstream(args.vcf)
    outstream = kevlar.open(args.out, 'w')
    writer = kevlar.vcf.VCFWriter(outstream, source='kevlar::varfilter')
    writer.write_header()
    varmask = load_variant_mask(kevlar.open(args.bed, 'r'))
    for varcall in varfilter(reader, varmask):
        writer.write(call)
