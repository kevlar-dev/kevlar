#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019 Battelle National Biodefense Institute
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar


def load_predictions(varcalls):
    message = 'Loading predictions to filter'
    kevlar.plog('[kevlar::varfilter]', message)
    progress_indictator = kevlar.ProgressIndicator(
        '[kevlar::varfilter]     {counter} intervals loaded', interval=1e4,
        breaks=[1e5, 1e6, 1e7], usetimer=True,
    )
    index = kevlar.IntervalForest()
    for call in varcalls:
        index.insert(*call.region, data=call)
        progress_indictator.update()
    return index


def varfilter(callstream, maskstream):
    callindex = load_predictions(callstream)
    message = 'Filtering preliminary variant calls'
    kevlar.plog('[kevlar::varfilter]', message)
    progress_indictator = kevlar.ProgressIndicator(
        '[kevlar::varfilter]     {counter} regions processed', interval=1e4,
        breaks=[1e5, 1e6, 1e7], usetimer=True,
    )
    for chrom, start, end, data in maskstream:
        hits = callindex.query(chrom, start, end)
        for interval in hits:
            interval.data.filter(kevlar.vcf.VariantFilter.UserFilter)
        progress_indictator.update()
    for varcall in callindex:
        yield varcall


def main(args):
    reader = kevlar.vcf.vcfstream(args.vcf)
    bedstream = kevlar.parse_bed(kevlar.open(args.filt, 'r'))
    outstream = kevlar.open(args.out, 'w')
    writer = kevlar.vcf.VCFWriter(outstream, source='kevlar::varfilter')
    writer.write_header()
    for varcall in varfilter(reader, bedstream):
        writer.write(varcall)
