#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
import re


def cutout(partstream, refrfile, **kwargs):
    for partition in partstream:
        contigs = list(partition)
        ccmatch = re.search(r'kvcc=(\d+)', contigs[0].name)
        cc = ccmatch.group(1) if ccmatch else None
        for gdna in kevlar.localize.localize(contigs, refrfile, **kwargs):
            yield cc, gdna


def main(args):
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(args.contigs, 'r'))
    if args.part_id:
        pstream = kevlar.parse_single_partition(contigstream, args.part_id)
    else:
        pstream = kevlar.parse_partitioned_reads(contigstream)
    outstream = kevlar.open(args.out, 'w')
    localizer = cutout(
        pstream, args.refr, seedsize=args.seed_size, delta=args.delta,
        maxdiff=args.max_diff, inclpattern=args.include,
        exclpattern=args.exclude, logstream=args.logfile
    )
    for part, gdna in localizer:
        seqname = gdna.defline
        if part is not None:
            seqname += ' kvcc={}'.format(part)
        record = kevlar.sequence.Record(name=seqname, sequence=gdna.sequence)
        kevlar.sequence.write_record(record, outstream)
