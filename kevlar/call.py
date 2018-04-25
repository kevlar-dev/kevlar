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


def align_mates(record, refrfile):
    fasta = ''
    for n, mateseq in enumerate(record.mateseqs, 1):
        fasta += '>mateseq{:d}\n{:s}\n'.format(n, mateseq)
    cmd = 'bwa mem {:s} -'.format(refrfile)
    cmdargs = cmd.split()
    for seqid, start, end in bwa_align(cmdargs, seqstring=fasta):
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


def prelim_call(targetlist, querylist, match=1, mismatch=2, gapopen=5,
                gapextend=0, ksize=31, refrfile=None, debug=False, mindist=5,
                logstream=sys.stderr):
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
            if refrfile and len(query.mateseqs) > 0:
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
                if alignment.matedist:
                    avgdistance = '{:.2f}'.format(alignment.matedist)
                    varcall.annotate('MATEDIST', avgdistance)
                yield varcall


def call(*args, **kwargs):
    """Thin wrapper over the `prelim_call` function.

    This function applies a deduplication procedure to preliminary calls.
    """
    for call in dedup(prelim_call(*args, **kwargs)):
        yield call


def main(args):
    outstream = kevlar.open(args.out, 'w')
    qinstream = kevlar.parse_augmented_fastx(kevlar.open(args.queryseq, 'r'))
    queryseqs = list(qinstream)
    tinstream = kevlar.open(args.targetseq, 'r')
    targetseqs = list(kevlar.reference.load_refr_cutouts(tinstream))
    caller = call(
        targetseqs, queryseqs,
        args.match, args.mismatch, args.open, args.extend,
        args.ksize, args.refr, args.debug, 5, args.logfile
    )
    for varcall in caller:
        print(varcall.vcf, file=outstream)
