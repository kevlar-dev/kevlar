#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
import sys
import kevlar
from kevlar.reference import bwa_align
from kevlar.varmap import VariantMapping
from kevlar.vcf import VariantFilter as vf
import khmer
from khmer import _buckets_per_byte


def align_mates(record, refrfile):
    fasta = ''
    for n, mateseq in enumerate(record.mates, 1):
        fasta += '>mateseq{:d}\n{:s}\n'.format(n, mateseq)
    cmd = 'bwa mem {:s} -'.format(refrfile)
    cmdargs = cmd.split()
    for seqid, start, end, seq in bwa_align(cmdargs, seqstring=fasta):
        yield seqid, start, end


def mate_distance(mate_positions, gdna_position):
    gdnaseq, startpos, endpos = gdna_position

    def intvldist(istart, iend):
        x, y = sorted(((startpos, endpos), (istart, iend)))
        if x[0] <= x[1] < y[0]:
            return y[0] - x[1]
        return 0

    distances = list()
    for seqid, start, end in mate_positions:
        if seqid != gdnaseq:
            continue
        d = intvldist(start, end)
        if d < 10000:
            distances.append(d)
    if len(distances) == 0:
        return float('Inf')
    return sum(distances) / len(distances)


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
    for seqid in calls:
        for position in calls[seqid]:
            sortedcalls = sorted(
                calls[seqid][position], key=lambda call: call.windowlength,
                reverse=True
            )
            yield sortedcalls[0]


def prelim_call(targetlist, querylist, partid=None, match=1, mismatch=2,
                gapopen=5, gapextend=0, ksize=31, refrfile=None, debug=False,
                mindist=5, logstream=sys.stderr):
    """Implement the `kevlar call` procedure as a generator function."""
    for query in sorted(querylist, reverse=True, key=len):
        alignments = list()
        for target in sorted(targetlist, key=lambda cutout: cutout.defline):
            mapping = VariantMapping(
                query, target, match=match, mismatch=mismatch, gapopen=gapopen,
                gapextend=gapextend
            )
            alignments.append(mapping)
        aligns2report = alignments_to_report(alignments)
        if len(aligns2report) > 1:
            if refrfile and len(query.mates) > 0:
                mate_pos = list(align_mates(query, refrfile))
                if len(mate_pos) > 0:
                    for aln in aligns2report:
                        aln.matedist = mate_distance(mate_pos, aln.interval)
                    aligns2report.sort(key=lambda aln: aln.matedist)

        for n, alignment in enumerate(aligns2report):
            if debug:
                print('DEBUG ', alignment.cutout.defline, ' vs ',
                      alignment.contig.name, '\n', str(alignment), sep='',
                      end='\n\n', file=logstream)
            for varcall in alignment.call_variants(ksize, mindist, logstream):
                if partid is not None:
                    varcall.annotate('PART', partid)
                if alignment.matedist:
                    varcall.annotate('MATEDIST', alignment.matedist)
                yield varcall


def call(*args, **kwargs):
    """Thin wrapper over the `prelim_call` function.

    This function applies a deduplication procedure to preliminary calls.
    """
    for call in dedup(prelim_call(*args, **kwargs)):
        yield call


def load_contigs(contigstream, logstream=sys.stderr):
    message = 'loading contigs into memory by partition'
    print('[kevlar::call]', message, file=logstream)
    contigs_by_partition = dict()
    nparts = 0
    ncontigs = 0
    for partid, contiglist in contigstream:
        nparts += 1
        ncontigs += len(contiglist)
        contigs_by_partition[partid] = contiglist
    message = 'loaded {} contigs from {} partitions'.format(ncontigs, nparts)
    print('[kevlar::call]', message, file=logstream)
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
        print('[kevlar::call]', message, file=args.logfile)
        ntables = 4
        buckets = args.mask_mem * _buckets_per_byte['nodegraph'] / ntables
        mask = khmer.Nodetable(args.ksize, buckets, ntables)
    progress_indicator = kevlar.ProgressIndicator(
        '[kevlar::call] processed contigs/gDNAs for {counter} partitions',
        interval=10, breaks=[100, 1000, 10000], logstream=args.logfile,
    )
    for partid, gdnas in gdnastream:
        progress_indicator.update()
        contigs = contigs_by_partition[partid]
        caller = call(
            gdnas, contigs, partid, match=args.match, mismatch=args.mismatch,
            gapopen=args.open, gapextend=args.extend, ksize=args.ksize,
            refrfile=args.refr, debug=args.debug, mindist=5,
            logstream=args.logfile
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
            print('[kevlar::call]', message, file=args.logfile)
        mask.save(args.gen_mask)
