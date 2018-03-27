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
from kevlar.varmap import VariantMapping


def align_mates(record, refrfile):
    fasta = ''
    for n, mateseq in enumerate(record.mateseqs, 1):
        fasta += '>mateseq{:d}\n{:s}\n'.format(n, mateseq)
    cmd = 'bwa mem {:s} -'.format(refrfile)
    cmdargs = cmd.split()
    for seqid, pos in kevlar.bwa_align(cmdargs, seqstring=fasta):
        yield seqid, pos


def mate_distance(mate_positions, gdna_position):
    gdnaseq, startpos, endpos = gdna_position

    def pointdist(point):
        if point < startpos:
            return startpos - point
        elif point > endpos:
            return point - endpos
        else:
            return 0

    distances = list()
    for seqid, pos in mate_positions:
        if seqid != gdnaseq:
            continue
        d = pointdist(pos)
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


def call(targetlist, querylist, match=1, mismatch=2, gapopen=5,
         gapextend=0, ksize=31, refrfile=None, mindist=5,
         logstream=sys.stderr):
    """Wrap the `kevlar call` procedure as a generator function.

    Input is the following.
    - an iterable containing one or more target sequences from the reference
      genome, stored as khmer or screed sequence records
    - an iterable containing one or more contigs assembled by kevlar, stored as
      khmer or screed sequence records
    - alignment match score (integer)
    - alignment mismatch penalty (integer)
    - alignment gap open penalty (integer)
    - alignment gap extension penalty (integer)
    - mates of interesting reads, in case these are needed to distinguish
      between multiple best hist (filename)
    - reference file to which mates of interesting reads, if any, will be
      mapped to disambiguate multi-mapping contigs

    The function yields tuples of target sequence name, query sequence name,
    and alignment CIGAR string
    """
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
            for varcall in alignment.call_variants(ksize, mindist, logstream):
                if alignment.matedist:
                    varcall.info['MD'] = '{:.2f}'.format(alignment.matedist)
                    if n > 0:
                        varcall.annotate('NC', 'matefail')
                yield varcall


def main(args):
    outstream = kevlar.open(args.out, 'w')
    qinstream = kevlar.parse_augmented_fastx(kevlar.open(args.queryseq, 'r'))
    queryseqs = list(qinstream)
    tinstream = kevlar.open(args.targetseq, 'r')
    targetseqs = list(kevlar.reference.load_refr_cutouts(tinstream))
    caller = call(
        targetseqs, queryseqs,
        args.match, args.mismatch, args.open, args.extend,
        args.ksize, args.refr
    )
    for varcall in caller:
        print(varcall.vcf, file=outstream)
