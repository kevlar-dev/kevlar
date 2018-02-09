#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import kevlar
from kevlar.assemble import assemble_greedy, assemble_fml_asm
from kevlar.augment import augment
from kevlar.localize import localize
from kevlar.call import call


def augment_and_mark(augseqs, nakedseqs):
    for seq in augment(augseqs, nakedseqs):
        seq.name = '{:s} nikmers={:d}'.format(seq.name, len(seq.ikmers))
        yield seq


def alac(pstream, refrfile, ksize=31, bigpart=10000, delta=25, seedsize=31,
         maxdiff=10000, match=1, mismatch=2, gapopen=5, gapextend=0,
         greedy=False, min_ikmers=None, logstream=sys.stderr):
    assembler = assemble_greedy if greedy else assemble_fml_asm
    for partition in pstream:
        reads = list(partition)
        if len(reads) > bigpart:
            message = 'skipping partition with {:d} reads'.format(len(reads))
            print('[kevlar::alac] WARNING:', message, file=logstream)
            continue

        # Assemble partitioned reads into contig(s)
        contigs = list(assembler(reads, logstream=logstream))
        contigs = list(augment_and_mark(reads, contigs))
        if min_ikmers is not None:
            # Apply min ikmer filter if it's set
            contigs = [c for c in contigs if len(c.ikmers) >= min_ikmers]
        if len(contigs) == 0:
            continue

        # Identify the genomic region(s) associated with each contig
        localizer = localize(contigs, refrfile, seedsize, delta=delta,
                             logstream=logstream)
        targets = list(localizer)
        if len(targets) == 0:
            continue

        # Align contigs to genomic targets to make variant calls
        caller = call(targets, contigs, match, mismatch, gapopen, gapextend,
                      ksize)
        for varcall in caller:
            yield varcall


def main(args):
    readstream = kevlar.parse_augmented_fastx(kevlar.open(args.infile, 'r'))
    if args.part_id:
        pstream = kevlar.parse_single_partition(readstream, args.part_id)
    else:
        pstream = kevlar.parse_partitioned_reads(readstream)
    outstream = kevlar.open(args.out, 'w')
    workflow = alac(
        pstream, args.refr, ksize=args.ksize, bigpart=args.bigpart,
        delta=args.delta, seedsize=args.seed_size, maxdiff=args.max_diff,
        match=args.match, mismatch=args.mismatch, gapopen=args.open,
        gapextend=args.extend, greedy=args.greedy, min_ikmers=args.min_ikmers,
        logstream=args.logfile
    )

    for varcall in workflow:
        print(varcall.vcf, file=outstream)
