#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import re
import sys
import kevlar
from kevlar.assemble import assemble_greedy, assemble_fml_asm
from kevlar.localize import localize
from kevlar.call import call


def alac(pstream, refrfile, ksize=31, bigpart=10000, delta=50, seedsize=31,
         maxdiff=None, match=1, mismatch=2, gapopen=5, gapextend=0,
         greedy=False, fallback=False, min_ikmers=None, logstream=sys.stderr):
    assembler = assemble_greedy if greedy else assemble_fml_asm
    for partition in pstream:
        reads = list(partition)
        ccmatch = re.search(r'kvcc=(\d+)', reads[0].name)
        if len(reads) > bigpart:
            message = 'skipping partition with {:d} reads'.format(len(reads))
            print('[kevlar::alac] WARNING:', message, file=logstream)
            continue

        # Assemble partitioned reads into contig(s)
        contigs = list(assembler(reads, logstream=logstream))
        if len(contigs) == 0 and assembler == assemble_fml_asm and fallback:
            message = 'WARNING: no contig assembled by fermi-lite'
            if ccmatch:
                message += ' for CC={:s}'.format(ccmatch.group(1))
            message += '; attempting again with home-grown greedy assembler'
            print('[kevlar::alac]', message, file=logstream)
            contigs = list(assemble_greedy(reads, logstream=logstream))
        if min_ikmers is not None:
            # Apply min ikmer filter if it's set
            contigs = [c for c in contigs if len(c.ikmers) >= min_ikmers]
        if len(contigs) == 0:
            continue

        # Identify the genomic region(s) associated with each contig
        refrstream = kevlar.open(refrfile, 'r')
        seqs = kevlar.seqio.parse_seq_dict(refrstream)
        localizer = localize(contigs, refrfile, seedsize, delta=delta,
                             maxdiff=maxdiff, refrseqs=seqs,
                             logstream=logstream)
        targets = list(localizer)
        if len(targets) == 0:
            continue

        # Align contigs to genomic targets to make variant calls
        caller = call(targets, contigs, match, mismatch, gapopen, gapextend,
                      ksize, refrfile)
        for varcall in caller:
            if ccmatch:
                varcall.annotate('PART', ccmatch.group(1))
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
        gapextend=args.extend, greedy=args.greedy, fallback=args.fallback,
        min_ikmers=args.min_ikmers, logstream=args.logfile
    )

    for varcall in workflow:
        print(varcall.vcf, file=outstream)
