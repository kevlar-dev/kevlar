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
from kevlar.reference import bwa_align
from kevlar.varmap import VariantMapping
from kevlar.vcf import VariantFilter as vf
import khmer
from khmer import _buckets_per_byte


def alignments_to_report(alignments):
    """Determine which alignments should be reported and used to call variants.

    In the simplest and best case, there is only a single alignment to
    consider. If there is more than one alignment, determine which ones are
    interpretable as a variant, and of these return the alignment(s) with the
    optimal score.
    """
    if len(alignments) == 1:
        return alignments
    scrtbl = [aln for aln in alignments if aln.vartype is not None]
    if len(scrtbl) == 0:
        finallist = alignments
    else:
        finallist = scrtbl
    bestscore = max([aln.score for aln in finallist])
    aligns2report = [aln for aln in finallist if aln.score == bestscore]
    return aligns2report


def dedup(callstream):
    calls = dict()
    for call in callstream:
        if call.seqid not in calls:
            calls[call.seqid] = defaultdict(set)
        calls[call.seqid][call.position].add(call)
    for seqid in sorted(calls):
        for position in sorted(calls[seqid]):
            sortedcalls = sorted(
                calls[seqid][position], key=lambda call: call.windowlength,
                reverse=True
            )
            yield sortedcalls[0]


def merge_adjacent(callstream):
    prev = None
    for call in callstream:
        if prev is not None:
            trymerge = prev.test_merge(call)
            if trymerge is not None:
                call = trymerge
                prev = None
        if prev is not None:
            yield prev
        prev = call
    yield prev


def prelim_call(targetlist, querylist, partid=None, match=1, mismatch=2,
                gapopen=5, gapextend=0, ksize=31, refrfile=None, debug=False,
                mindist=5, homopolyfilt=True):
    """Implement the `kevlar call` procedure as a generator function."""
    for query in sorted(querylist, reverse=True, key=len):
        alignments = list()
        for target in sorted(targetlist, key=lambda cutout: cutout.defline):
            mapping = VariantMapping(
                query, target, match=match, mismatch=mismatch, gapopen=gapopen,
                gapextend=gapextend, homopolyfilt=homopolyfilt,
            )
            alignments.append(mapping)
        aligns2report = alignments_to_report(alignments)
        for n, alignment in enumerate(aligns2report):
            if debug:
                kevlar.plog(
                    'DEBUG ', alignment.cutout.defline, ' vs ',
                    alignment.contig.name, '\n', str(alignment), sep='',
                    end='\n\n',
                )
            for varcall in alignment.call_variants(ksize, mindist):
                if partid is not None:
                    varcall.annotate('PART', partid)
                yield varcall


def call(*args, **kwargs):
    """Thin wrapper over the `prelim_call` function.

    This function applies a deduplication procedure to preliminary calls.
    """
    for call in merge_adjacent(dedup(prelim_call(*args, **kwargs))):
        yield call


def load_contigs(contigstream):
    kevlar.plog('[kevlar::call] Loading contigs into memory by partition')
    contigs_by_partition = dict()
    nparts = 0
    ncontigs = 0
    for partid, contiglist in contigstream:
        nparts += 1
        ncontigs += len(contiglist)
        contigs_by_partition[partid] = contiglist
    message = 'Loaded {} contigs from {} partitions'.format(ncontigs, nparts)
    kevlar.plog('[kevlar::call]', message)
    return contigs_by_partition


def main(args):
    # Input and output files
    outstream = kevlar.open(args.out, 'w')
    writer = kevlar.vcf.VCFWriter(
        outstream, source='kevlar::call', refr=args.refr,
    )
    writer.write_header()

    # Contigs = query sequences
    contigstream = kevlar.parse_partitioned_reads(
        kevlar.parse_augmented_fastx(kevlar.open(args.queryseq, 'r'))
    )
    contigs_by_partition = load_contigs(contigstream)

    gdnastream = kevlar.parse_partitioned_reads(
        kevlar.reference.load_refr_cutouts(kevlar.open(args.targetseq, 'r'))
    )
    mask = None
    if args.gen_mask:
        message = 'generating mask of variant-spanning k-mers'
        kevlar.plog('[kevlar::call]', message)
        ntables = 4
        buckets = args.mask_mem * _buckets_per_byte['nodegraph'] / ntables
        mask = khmer.Nodetable(args.ksize, buckets, ntables)
    progress_indicator = kevlar.ProgressIndicator(
        '[kevlar::call] processed contigs/gDNAs for {counter} partitions',
        interval=10, breaks=[100, 1000, 10000],
    )
    for partid, gdnas in gdnastream:
        progress_indicator.update()
        if partid not in contigs_by_partition:
            continue
        contigs = contigs_by_partition[partid]
        caller = call(
            gdnas, contigs, partid, match=args.match, mismatch=args.mismatch,
            gapopen=args.open, gapextend=args.extend, ksize=args.ksize,
            refrfile=args.refr, debug=args.debug, mindist=5,
            homopolyfilt=not args.no_homopoly_filter,
        )
        for varcall in caller:
            if args.gen_mask:
                window = varcall.attribute('ALTWINDOW')
                if window is not None and len(window) >= args.ksize:
                    mask.consume(window)
            writer.write(varcall)
    if args.gen_mask:
        fpr = khmer.calc_expected_collisions(mask, max_false_pos=1.0)
        if fpr > args.mask_max_fpr:
            message = 'WARNING: mask FPR is {:.4f}'.format(fpr)
            message += '; exceeds user-specified limit'
            message += ' of {:.4f}'.format(args.mask_max_fpr)
            kevlar.plog('[kevlar::call]', message)
        mask.save(args.gen_mask)
