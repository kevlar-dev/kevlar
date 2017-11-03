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
from kevlar.novel import novel, load_samples
from kevlar.filter import filter as kfilter
from kevlar.filter import load_mask
from kevlar.partition import partition
from kevlar.alac import alac


def simplex(case, casecounts, controlcounts, refrfile, ctrlmax=0, casemin=5,
            mask=None, filtermem=1e6, filterfpr=0.001,
            partminabund=2, partmaxabund=200, dedup=True,
            delta=50, match=1, mismatch=2, gapopen=5, gapextend=0,
            ksize=31, logstream=sys.stderr):
    """
    Execute the simplex germline variant discovery workflow.

    Parameters for identifying novel k-mers:
    - case: stream of input reads from case sample
    - casecounts: a Counttable object containing k-mer counts from the case
                  sample
    - controlcounts: a list of Counttable objects containing k-mer counts from
                     the control samples
    - refrfile: BWA-indexed reference genome sequences
    - ctrlmax: maximum abundance in each control sample for a k-mer to be
               designated "interesting"
    - casemin: minimum abundance in the case sample for a k-mer to be
               designated "interesting"

    Parameters for filtering "interesting" reads:
    - mask: Nodetable containing k-mers from reference sequence, contaminants
    - filtermem: memory to allocate for recalculating k-mer abundances
    - filterfpr: abort if FDR of recomputed k-mer abundances is too high

    Parameters for partitioning reads into distinct sets
    - partminabund: discard an "interesting" k-mer if it is annotated in < this
                    many reads
    - partmaxabund: discard an "interesting" k-mer if it is annotated in > this
                    many reads
    - dedup: boolean indicating whether PCR duplicate removal should be run

    Parameters for assembling, aligning, and calling variants:
    - delta: number of bp to extend genomic cutout
    - match: alignment match score
    - mismatch: alignment mismatch penalty
    - gapopen: alignment gap open penalty
    - gapextend: alignment gap extension penalty
    """
    discoverer = novel(
        case, [casecounts], controlcounts, ksize=ksize, casemin=casemin,
        ctrlmax=ctrlmax, updateint=100000, logstream=logstream
    )
    filterer = kfilter(
        discoverer, mask=mask, minabund=casemin, ksize=ksize,
        memory=filtermem, maxfpr=filterfpr, logstream=logstream
    )
    partitioner = partition(
        filterer, dedup=dedup, minabund=partminabund, maxabund=partmaxabund,
        logstream=logstream
    )
    for cc in partitioner:
        caller = alac(
            cc, refrfile, ksize=ksize, delta=delta, match=match,
            mismatch=mismatch, gapopen=gapopen, gapextend=gapextend,
            logstream=logstream
        )
        for variant in caller:
            yield variant


def main(args):
    cases = load_samples(
        args.case_counts, [args.case], ksize=args.ksize,
        memory=args.novel_memory, maxfpr=args.novel_fpr,
        numthreads=args.threads, logstream=args.logfile
    )
    controls = load_samples(
        args.control_counts, args.control, ksize=args.ksize,
        memory=args.novel_memory, maxfpr=args.novel_fpr,
        numthreads=args.threads, logstream=args.logfile
    )
    mask = load_mask(
        args.mask_files, args.ksize, args.mask_memory, maxfpr=args.filter_fpr,
        logstream=args.logfile
    )

    outstream = kevlar.open(args.out, 'w')
    caserecords = kevlar.multi_file_iter_screed(args.case)
    workflow = simplex(
        caserecords, cases[0], controls, args.refr, ksize=args.ksize,
        ctrlmax=args.ctrl_max, casemin=args.case_min, mask=mask,
        filtermem=args.filter_memory, filterfpr=args.filter_fpr,
        partminabund=args.part_min_abund, partmaxabund=args.part_max_abund,
        dedup=args.dedup, delta=args.delta, match=args.match,
        mismatch=args.mismatch, gapopen=args.open, gapextend=args.extend,
        logstream=args.logfile
    )
    for variant in workflow:
        print(variant.vcf, file=outstream)
