#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import khmer
import kevlar
from kevlar.novel import novel, load_samples
from kevlar.filter import filter as kfilter
from kevlar.filter import load_mask
from kevlar.partition import partition
from kevlar.alac import alac
from kevlar.simlike import simlike


def populate_refrsct(refrseqfile, refrsctfile=None, refrsctmem=1e9, ksize=31,
                     logstream=sys.stderr):
    if refrsctfile:
        return khmer.SmallCounttable.load(refrsctfile)
    else:
        sct = khmer.SmallCounttable(ksize, int(refrsctmem) / 4, 4)
        sct.consume_seqfile(refrseqfile)
        fpr = kevlar.sketch.estimate_fpr(sct)
        if fpr > 0.05:
            message = 'WARNING: false positive rate of {:.4f}'.format(fpr)
            message += ' may be problematic for likelihood calculations.'
            message += ' Consider increasing memory for reference '
            message += ' k-mer counts.'
            print('[kevlar::simplex]', warning, file=logstream)
        return sct


def simplex(case, casecounts, controlcounts, refrfile, ctrlmax=0, casemin=5,
            mask=None, filtermem=1e6, filterfpr=0.001,
            partminabund=2, partmaxabund=200, dedup=True,
            delta=50, seedsize=31, match=1, mismatch=2, gapopen=5, gapextend=0,
            refrsctfile=None, refrsctmem=1e9, mu=30.0,
            sigma=8.0, epsilon=0.001, labels=None, ksize=31, threads=1,
            logstream=sys.stderr):
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
    - seedsize: size of seeds to use for localizing contigs
    - delta: number of bp to extend genomic cutout
    - match: alignment match score
    - mismatch: alignment mismatch penalty
    - gapopen: alignment gap open penalty
    - gapextend: alignment gap extension penalty

    Parameters for computing likelihood scores
    - refrsctfile: smallcounttable of k-mer abundances in the reference genome;
                   if not provided, will be populated from the refrfile
                   parameter
    - refrsctmem: memory to allocate for reference smallcounttable (if it needs
                  to be populated from scratch)
    - mu: observed/expected average k-mer coverage
    - sigma: observed/expected standard deviation for k-mer coverage
    - epsilon: base error rate
    - labels: human-readable lables for each sample, case/proband first
    """
    discoverer = novel(
        case, [casecounts], controlcounts, ksize=ksize, casemin=casemin,
        ctrlmax=ctrlmax, logstream=logstream
    )
    filterer = kfilter(
        discoverer, mask=mask, casemin=casemin, ctrlmax=ctrlmax, ksize=ksize,
        memory=filtermem, maxfpr=filterfpr, logstream=logstream
    )
    partitioner = partition(
        filterer, dedup=dedup, minabund=partminabund, maxabund=partmaxabund,
        logstream=logstream
    )

    caller = alac(
        partitioner, refrfile, threads=threads, ksize=ksize, delta=delta,
        seedsize=seedsize, match=match, mismatch=mismatch, gapopen=gapopen,
        gapextend=gapextend, logstream=logstream
    )

    refrsct = populate_refrsct(refrfile, refrsctfile, refrsctmem=refrsctmem,
                               ksize=ksize, logstream=logstream)
    scorer = simlike(
        caller, casecounts, controlcounts, refrsct, mu=mu, sigma=sigma,
        epsilon=epsilon, casemin=casemin, samplelabels=labels,
        logstream=logstream
    )

    for variant in scorer:
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
    kevlar.novel.save_all_counts(
        cases, args.save_case_counts, controls, args.save_ctrl_counts,
        logstream=args.logfile
    )

    if not args.labels:
        args.labels = kevlar.simlike.default_sample_labels(len(controls) + 1)
    outstream = kevlar.open(args.out, 'w')
    caserecords = kevlar.multi_file_iter_screed(args.case)
    workflow = simplex(
        caserecords, cases[0], controls, args.refr, ksize=args.ksize,
        ctrlmax=args.ctrl_max, casemin=args.case_min, mask=mask,
        filtermem=args.filter_memory, filterfpr=args.filter_fpr,
        partminabund=args.part_min_abund, partmaxabund=args.part_max_abund,
        dedup=args.dedup, delta=args.delta, seedsize=args.seed_size,
        match=args.match, mismatch=args.mismatch, gapopen=args.open,
        gapextend=args.extend, threads=args.threads, refrsctfile=args.refr_sct,
        refrsctmem=args.refr_sct_mem, mu=args.mu, sigma=args.sigma,
        epsilon=args.epsilon, labels=args.labels, logstream=args.logfile
    )
    writer = kevlar.vcf.VCFWriter(
        outstream, source='kevlar::simplex', refr=args.refr,
    )
    for label in args.labels:
        writer.register_sample(label)
    writer.write_header()
    for varcall in workflow:
        writer.write(varcall)
