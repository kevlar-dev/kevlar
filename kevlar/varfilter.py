#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019 Battelle National Biodefense Institute
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
import pickle


def save_index_to_file(index, filename):
    with open(filename, 'wb') as fh:
        pickle.dump(index, fh)


def load_index_from_file(filename):
    with open(filename, 'rb') as fh:
        index = pickle.load(fh)
    return index


def load_variant_mask(bedstream):
    message = 'Loading genomic intervals to be filtered into memory'
    kevlar.plog('[kevlar::varfilter]', message)
    progress_indictator = kevlar.ProgressIndicator(
        '[kevlar::varfilter]     {counter} intervals loaded', interval=1e4,
        breaks=[1e5, 1e6, 1e7], usetimer=True,
    )
    index = kevlar.IntervalForest()
    for chrom, start, end, data in kevlar.parse_bed(bedstream):
        index.insert(chrom, start, end)
        progress_indictator.update()
    return index


def varfilter(calls, varmask):
    message = 'Filtering preliminary variant calls'
    kevlar.plog('[kevlar::varfilter]', message)
    progress_indictator = kevlar.ProgressIndicator(
        '[kevlar::varfilter]     {counter} calls processed', interval=1e3,
        breaks=[1e4, 1e5, 1e6], usetimer=True,
    )
    for varcall in calls:
        if varmask.query(varcall.seqid, varcall.position) != set():
            varcall.filter(kevlar.vcf.VariantFilter.UserFilter)
        yield varcall
        progress_indictator.update()


def main(args):
    reader = kevlar.vcf.vcfstream(args.vcf)
    outstream = kevlar.open(args.out, 'w')
    writer = kevlar.vcf.VCFWriter(outstream, source='kevlar::varfilter')
    writer.write_header()
    if args.load_index:
        varmask = load_index_from_file(args.filt)
    else:
        varmask = load_variant_mask(kevlar.open(args.filt, 'r'))
    if args.save_index:
        save_index_to_file(varmask, args.save_index)
    for varcall in varfilter(reader, varmask):
        writer.write(varcall)
