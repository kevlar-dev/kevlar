#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from math import log
import re
import sys
import khmer
import kevlar
from kevlar.vcf import Variant


class VariantMapping(object):
    def __init__(self, contig, cutout, score, cigar, strand=1):
        self.contig = contig
        self.cutout = cutout
        self.score = score
        self.cigar = cigar
        self.strand = strand
        self.matedist = None

    @property
    def interval(self):
        return self.cutout.interval

    @property
    def varseq(self):
        assert self.strand in (-1, 1)
        if self.strand == 1:
            return self.contig.sequence
        else:
            return kevlar.revcom(self.contig.sequence)

    @property
    def refrseq(self):
        return self.cutout.sequence

    @property
    def seqid(self):
        return self.cutout._seqid

    @property
    def pos(self):
        return self.cutout._startpos

    def call_variants(self, ksize):
        snvmatch = re.search('^(\d+)([DI])(\d+)M(\d+)[DI]$', self.cigar)
        snvmatch2 = re.search('^(\d+)([DI])(\d+)M(\d+)[DI](\d+)M$', self.cigar)
        if snvmatch:
            offset = int(snvmatch.group(1))
            if snvmatch.group(2) == 'I':
                offset *= -1
            length = int(snvmatch.group(3))
            return call_snv(self, offset, length, ksize)
        elif snvmatch2 and int(snvmatch2.group(5)) <= 5:
            offset = int(snvmatch2.group(1))
            if snvmatch2.group(2) == 'I':
                offset *= -1
            length = int(snvmatch2.group(3))
            return call_snv(self, offset, length, ksize)

        indelmatch = re.search(
            '^(\d+)([DI])(\d+)M(\d+)([ID])(\d+)M(\d+)[DI]$', self.cigar
        )
        indelmatch2 = re.search(
            '^(\d+)([DI])(\d+)M(\d+)([ID])(\d+)M(\d+)[DI](\d+)M$', self.cigar
        )
        if indelmatch:
            offset = int(indelmatch.group(1))
            if indelmatch.group(2) == 'I':
                offset *= -1
            leftmatch = int(indelmatch.group(3))
            indellength = int(indelmatch.group(4))
            indeltype = indelmatch.group(5)
            callfunc = call_deletion if indeltype == 'D' else call_insertion
            return callfunc(self, offset, ksize, leftmatch, indellength)
        elif indelmatch2 and int(indelmatch2.group(8)) <= 5:
            offset = int(indelmatch2.group(1))
            if indelmatch2.group(2) == 'I':
                offset *= -1
            leftmatch = int(indelmatch2.group(3))
            indellength = int(indelmatch2.group(4))
            indeltype = indelmatch2.group(5)
            callfunc = call_deletion if indeltype == 'D' else call_insertion
            return callfunc(self, offset, ksize, leftmatch, indellength)

        nocall = Variant(
            self.seqid, self.pos, '.', '.', NC='inscrutablecigar',
            CS=self.contig.sequence, CIGAR=self.cigar,
        )
        return [nocall]


def call_snv(aln, offset, length, ksize):
    targetshort = False
    if offset < 0:
        offset *= -1
        gdnaoffset = 0
        targetshort = True
        t = aln.refrseq[:length]
        q = aln.varseq[offset:offset+length]
    else:
        gdnaoffset = offset
        t = aln.refrseq[offset:offset+length]
        q = aln.varseq[:length]
    diffs = [(i, t[i], q[i]) for i in range(length) if t[i] != q[i]]
    if len(diffs) == 0:
        nocall = Variant(aln.seqid, aln.cutout.local_to_global(gdnaoffset),
                         '.', '.', NC='perfectmatch', QN=aln.contig.name, QS=q)
        return [nocall]

    snvs = list()
    for diff in diffs:
        minpos = max(diff[0] - ksize + 1, 0)
        maxpos = min(diff[0] + ksize, length)
        window = q[minpos:maxpos]
        refrwindow = t[minpos:maxpos]

        # numoverlappingkmers = len(window) - ksize + 1
        # kmers = [window[i:i+ksize] for i in range(numoverlappingkmers)]
        refr = diff[1].upper()
        alt = diff[2].upper()
        localcoord = diff[0]
        if not targetshort:
            localcoord += offset
        globalcoord = aln.cutout.local_to_global(localcoord)
        snv = Variant(aln.seqid, globalcoord, refr, alt, VW=window,
                      RW=refrwindow, IK=str(len(aln.contig.ikmers)))
        snvs.append(snv)
    return snvs


def deletion_allele(target, query, offset, ksize, leftmatch, indellength):
    minpos = leftmatch - ksize + 1
    maxpos = leftmatch + ksize - 1
    altwindow = query[minpos:maxpos]
    minpos += offset
    maxpos += offset + indellength
    refrwindow = target[minpos:maxpos]

    refr = target[offset+leftmatch-1:offset+leftmatch+indellength]
    alt = refr[0]
    return refr, alt, refrwindow, altwindow


def insertion_allele(target, query, offset, ksize, leftmatch, indellength):
    minpos = leftmatch - ksize + 1
    maxpos = leftmatch + ksize + indellength - 1
    altwindow = query[minpos:maxpos]
    minpos += offset
    maxpos += offset - indellength
    refrwindow = target[minpos:maxpos]

    alt = query[leftmatch-1:leftmatch+indellength]
    refr = alt[0]
    return refr, alt, refrwindow, altwindow


def call_deletion(aln, offset, ksize, leftmatch, indellength):
    if offset < 0:
        offset *= -1
        targetshort = True
        alt, refr, altwindow, refrwindow = insertion_allele(
            aln.varseq, aln.refrseq, offset, ksize, leftmatch, indellength
        )
    else:
        targetshort = False
        refr, alt, refrwindow, altwindow = deletion_allele(
            aln.refrseq, aln.varseq, offset, ksize, leftmatch, indellength
        )
    # This assertion is no longer valid when query is longer than target
    # assert len(refr) == indellength + 1
    localcoord = leftmatch
    if not targetshort:
        localcoord += offset
    globalcoord = aln.cutout.local_to_global(localcoord)
    indel = Variant(aln.seqid, globalcoord - 1, refr, alt, VW=altwindow,
                    RW=refrwindow, IK=str(len(aln.contig.ikmers)))
    return [indel]


def call_insertion(aln, offset, ksize, leftmatch, indellength):
    if offset < 0:
        offset *= -1
        targetshort = True
        alt, refr, altwindow, refrwindow = deletion_allele(
            aln.varseq, aln.refrseq, offset, ksize, leftmatch, indellength
        )
    else:
        targetshort = False
        refr, alt, refrwindow, altwindow = insertion_allele(
            aln.refrseq, aln.varseq, offset, ksize, leftmatch, indellength
        )

    # This assertion is no longer valid when query is longer than target
    # assert len(alt) == indellength + 1
    localcoord = leftmatch
    if not targetshort:
        localcoord += offset
    globalcoord = aln.cutout.local_to_global(localcoord)
    indel = Variant(aln.seqid, globalcoord - 1, refr, alt, VW=altwindow,
                    RW=refrwindow, IK=str(len(aln.contig.ikmers)))
    return [indel]


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


def align_both_strands(target, query, match=1, mismatch=2, gapopen=5,
                       gapextend=0):
    cigar1, score1 = kevlar.align(
        target.sequence, query.sequence, match, mismatch, gapopen, gapextend
    )
    cigar2, score2 = kevlar.align(
        target.sequence, kevlar.revcom(query.sequence), match, mismatch,
        gapopen, gapextend
    )
    if score2 > score1:
        cigar = cigar2
        score = score2
        strand = -1
    else:
        cigar = cigar1
        score = score1
        strand = 1
    return VariantMapping(query, target, score, cigar, strand)


def alignment_interpretable(cigar):
    patterns = [
        '^(\d+)([DI])(\d+)M(\d+)[DI]$',
        '^(\d+)([DI])(\d+)M(\d+)[DI](\d+)M$',
        '^(\d+)([DI])(\d+)M(\d+)([ID])(\d+)M(\d+)[DI]$',
        '^(\d+)([DI])(\d+)M(\d+)([ID])(\d+)M(\d+)[DI](\d+)M$',
    ]
    for pattern in patterns:
        if re.search(pattern, cigar) is not None:
            return True
    return False


def alignments_to_report(alignments):
    """Determine which alignments should be reported and used to call variants.

    In the simplest and best case, there is only a single alignment to
    consider. If there is more than one alignment, determine which ones are
    interpretable as a variant, and of these return the alignment(s) with the
    optimal score.
    """
    if len(alignments) == 1:
        return alignments
    scrtbl = [aln for aln in alignments if alignment_interpretable(aln.cigar)]
    if len(scrtbl) == 0:
        finallist = alignments
    else:
        finallist = scrtbl
    bestscore = max([aln.score for aln in finallist])
    aligns2report = [aln for aln in finallist if aln.score == bestscore]
    return aligns2report


def call(targetlist, querylist, match=1, mismatch=2, gapopen=5,
         gapextend=0, ksize=31, refrfile=None):
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
    #dolike = casecounts is not None and controlcounts is not None
    varcalls = list()
    for query in sorted(querylist, reverse=True, key=len):
        alignments = list()
        for target in sorted(targetlist, key=lambda cutout: cutout.defline):
            mapping = align_both_strands(target, query, match, mismatch,
                                         gapopen, gapextend)
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
            for varcall in alignment.call_variants(ksize):
                if alignment.matedist:
                    varcall.info['MD'] = '{:.2f}'.format(alignment.matedist)
                    if n > 0:
                        varcall.annotate('NC', 'matefail')
                yield varcall

    # if dolike:
    #     varcalls.sort(key=lambda v: float(v.info['DN']), reverse=True)
    # for v in varcalls:
    #     yield v


def main(args):
    outstream = kevlar.open(args.out, 'w')
    writer = kevlar.vcf.VCFWriter(source='kevlar::call')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(args.queryseq, 'r'))
    queryseqs = list(qinstream)
    tinstream = kevlar.open(args.targetseq, 'r')
    targetseqs = list(kevlar.reference.load_refr_cutouts(tinstream))
    caller = call(
        targetseqs, queryseqs,
        match=args.match, mismatch=args.mismatch, gapopen=args.open,
        gapextend=args.extend, ksize=args.ksize, refrfile=args.refridx
    )

    for varcall in caller:
        writer.write(varcall, fh=outstream)
