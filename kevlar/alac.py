#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
import kevlar
from kevlar.assemble import assemble_fml_asm
from kevlar.localize import localize
from kevlar.call import call
import khmer
import re
import sys


def alac(pstream, refrfile, threads=1, ksize=31, maxreads=10000, delta=50,
         seedsize=31, maxdiff=None, inclpattern=None, exclpattern=None,
         match=1, mismatch=2, gapopen=5, gapextend=0, min_ikmers=None,
         maskfile=None, maskmem=1e6, maskmaxfpr=0.01, logstream=sys.stderr):
    assembler = kevlar.assemble.assemble(pstream, maxreads=maxreads)
    contigs_by_partition = defaultdict(list)
    for partid, contig in assembler:
        contigs_by_partition[partid].append(contig)

    contigstream = [(pid, ctgs) for pid, ctgs in contigs_by_partition.items()]
    targeter = kevlar.localize.localize(
        contigstream, refrfile, seedsize=seedsize, delta=delta,
        maxdiff=maxdiff, inclpattern=inclpattern, exclpattern=exclpattern,
        logstream=logstream
    )
    targets_by_partition = defaultdict(list)
    for partid, gdna in targeter:
        targets_by_partition[partid].append(gdna)

    calls = list()
    for partid in sorted(targets_by_partition):
        gdnalist = targets_by_partition[partid]
        contigs = contigs_by_partition[partid]
        caller = kevlar.call.call(
            gdnalist, contigs, partid, match=match, mismatch=mismatch,
            gapopen=gapopen, gapextend=gapextend, ksize=ksize,
            refrfile=refrfile,
        )
        partcalls = list(caller)
        calls.extend(partcalls)
    calls = sorted(calls, key=lambda c: (c.seqid, c.position))
    if maskfile:
        message = 'generating mask of variant-spanning k-mers'
        print('[kevlar::alac]', message, file=logstream)
        numtables = 4
        buckets = maskmem * khmer._buckets_per_byte['nodegraph'] / numtables
        mask = khmer.Nodetable(ksize, buckets, numtables)
        for call in calls:
            window = call.attribute('ALTWINDOW')
            if window is not None and len(window) >= ksize:
                mask.consume(window)
        fpr = khmer.calc_expected_collisions(mask, max_false_pos=1.0)
        if fpr > maskmaxfpr:
            message = 'WARNING: mask FPR is {:.4f}'.format(fpr)
            message += '; exceeds user-specified limit'
            message += ' of {:.4f}'.format(maskmaxfpr)
            print('[kevlar::alac]', message, file=logstream)
        mask.save(maskfile)
    for call in calls:
        yield call


def main(args):
    readstream = kevlar.parse_augmented_fastx(kevlar.open(args.infile, 'r'))
    if args.part_id:
        pstream = kevlar.parse_single_partition(readstream, args.part_id)
    else:
        pstream = kevlar.parse_partitioned_reads(readstream)
    outstream = kevlar.open(args.out, 'w')
    workflow = alac(
        pstream, args.refr, threads=args.threads, ksize=args.ksize,
        maxreads=args.max_reads, delta=args.delta, seedsize=args.seed_size,
        maxdiff=args.max_diff, inclpattern=args.include,
        exclpattern=args.exclude, match=args.match, mismatch=args.mismatch,
        gapopen=args.open, gapextend=args.extend, min_ikmers=args.min_ikmers,
        maskfile=args.gen_mask, maskmem=args.mask_mem,
        maskmaxfpr=args.mask_max_fpr, logstream=args.logfile,
    )

    writer = kevlar.vcf.VCFWriter(
        outstream, source='kevlar::alac', refr=args.refr,
    )
    writer.write_header()
    for varcall in workflow:
        writer.write(varcall)
