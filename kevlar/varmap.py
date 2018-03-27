#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
from kevlar.alignment import align_both_strands
from kevlar.vcf import Variant, VariantSNV, VariantIndel
import re
import sys


patterns = {
    '^(\d+)([DI])(\d+)M(\d+)[DI]$': 'snv',
    '^(\d+)([DI])(\d+)M(\d+)[DI](\d+)M$': 'snv',
    '^(\d+)([DI])(\d+)M(\d+)([ID])(\d+)M(\d+)[DI]$': 'indel',
    '^(\d+)([DI])(\d+)M(\d+)([ID])(\d+)M(\d+)[DI](\d+)M$': 'indel',
}


class VariantMapping(object):
    """Class for managing contig alignments to reference genome.

    All variant calls are made from alignments of contigs (assembled from reads
    containing novel sequence content) to a cutout of the reference genome.
    This class manages relevant data and implements the variant calling
    procedures.
    """
    def __init__(self, contig, cutout, score=None, cigar=None, strand=1,
                 match=1, mismatch=2, gapopen=5, gapextend=0):
        if score is None:
            score, cigar, strand = align_both_strands(
                cutout, contig, match, mismatch, gapopen, gapextend
            )
        self.contig = contig
        self.cutout = cutout
        self.score = score
        self.cigar = cigar
        self.strand = strand
        self.matedist = None
        self.vartype = None
        self.alnmatch = None
        for pattern, vartype in patterns.items():
            matchobj = re.match(pattern, cigar)
            if matchobj:
                self.alnmatch = matchobj
                self.vartype = vartype
                break

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

    @property
    def offset(self):
        if self.alnmatch is None:
            return None
        return int(self.alnmatch.group(1))

    @property
    def targetshort(self):
        if self.alnmatch is None:
            return None
        return self.alnmatch.group(2) == 'I'

    @property
    def leftmatchlen(self):
        if self.alnmatch is None or self.vartype != 'indel':
            return None
        return int(self.alnmatch.group(3))

    @property
    def indellength(self):
        if self.alnmatch is None or self.vartype != 'indel':
            return None
        return int(self.alnmatch.group(4))

    @property
    def indeltype(self):
        if self.alnmatch is None or self.vartype != 'indel':
            return None
        return self.alnmatch.group(5)

    @property
    def rightmatchlen(self):
        if self.alnmatch is None or self.vartype != 'indel':
            return None
        return int(self.alnmatch.group(6))

    def call_variants(self, ksize, mindist=5, logstream=sys.stderr):
        """Attempt to call variants from this contig alignment.

        If the alignment CIGAR matches a known pattern, the appropriate caller
        is invoked (SNV or INDEL caller). If not, a "no call" is reported.
        """
        if self.vartype == 'snv':
            for call in self.call_snv(ksize, mindist, logstream=logstream):
                yield call
        elif self.vartype == 'indel':
            caller = self.call_insertion
            if self.indeltype == 'D':
                caller = self.call_deletion
            for call in caller(ksize):
                yield call
            for call in self.call_indel_snvs(ksize):
                yield call
        else:
            nocall = Variant(
                self.seqid, self.pos, '.', '.', NC='inscrutablecigar',
                QN=self.contig.name, QS=self.varseq, CG=self.cigar,
            )
            yield nocall

    def snv_variant(self, qseq, tseq, mismatches, offset, ksize):
        """Call SNVs from the aligned mismatched sequences.

        The `qseq` and `tseq` are strings containing query and target sequences
        of identical length; `mismatches` is a list of positions where `qseq`
        and `tseq` do not match; `offset` is the number of 5' nucleotides in
        the target not aligned to the query; and `ksize` is used to compute a
        window that spans all reference allele k-mers in `tseq` and all
        alternate allele k-mers in `qseq`.
        """
        length = len(qseq)
        assert len(tseq) == length
        for pos in mismatches:
            minpos = max(pos - ksize + 1, 0)
            maxpos = min(pos + ksize, length)
            altwindow = qseq[minpos:maxpos]
            refrwindow = tseq[minpos:maxpos]

            refr = tseq[pos].upper()
            alt = qseq[pos].upper()
            localcoord = pos + offset
            globalcoord = self.cutout.local_to_global(localcoord)
            nikmers = n_ikmers_present(self.contig.ikmers, altwindow)
            snv = VariantSNV(
                self.seqid, globalcoord, refr, alt, VW=altwindow,
                RW=refrwindow, IK=str(nikmers)
            )
            yield snv

    def call_snv(self, ksize, mindist=5, logstream=sys.stderr):
        """Call SNVs from the given alignment."""
        length = int(self.alnmatch.group(3))
        offset = self.offset
        if self.targetshort:
            gdnaoffset = 0
            t = self.refrseq[:length]
            q = self.varseq[offset:offset+length]
        else:
            gdnaoffset = offset
            t = self.refrseq[offset:offset+length]
            q = self.varseq[:length]
        diffs = [i for i in range(length) if t[i] != q[i]]
        if mindist:
            diffs = trim_terminal_snvs(diffs, length, mindist, logstream)
        if len(diffs) == 0:
            nocall = Variant(
                self.seqid, self.cutout.local_to_global(gdnaoffset), '.', '.',
                NC='perfectmatch', QN=self.contig.name, QS=q
            )
            yield nocall
            return

        for call in self.snv_variant(q, t, diffs, offset, ksize):
            yield call

    def deletion_allele(self, target, query, ksize):
        offset = self.offset
        leftmatch = self.leftmatchlen
        indellength = self.indellength

        minpos = max(leftmatch - ksize + 1, 0)
        maxpos = min(leftmatch + ksize - 1, len(query))
        altwindow = query[minpos:maxpos]
        minpos += offset
        maxpos += offset + indellength
        refrwindow = target[minpos:maxpos]

        refr = target[offset+leftmatch-1:offset+leftmatch+indellength]
        alt = refr[0]
        return refr, alt, refrwindow, altwindow

    def insertion_allele(self, target, query, ksize):
        offset = self.offset
        leftmatch = self.leftmatchlen
        indellength = self.indellength

        minpos = max(leftmatch - ksize + 1, 0)
        maxpos = min(leftmatch + ksize + indellength - 1, len(query))
        altwindow = query[minpos:maxpos]
        minpos += offset
        maxpos += offset - indellength
        refrwindow = target[minpos:maxpos]

        alt = query[leftmatch-1:leftmatch+indellength]
        refr = alt[0]
        return refr, alt, refrwindow, altwindow

    def call_deletion(self, ksize):
        if self.targetshort:
            alt, refr, altwindow, refrwindow = self.insertion_allele(
                self.varseq, self.refrseq, ksize
            )
        else:
            refr, alt, refrwindow, altwindow = self.deletion_allele(
                self.refrseq, self.varseq, ksize
            )
        # This assertion is no longer valid when query is longer than target
        # assert len(refr) == indellength + 1
        localcoord = self.leftmatchlen
        if not self.targetshort:
            localcoord += self.offset
        globalcoord = self.cutout.local_to_global(localcoord)
        nikmers = n_ikmers_present(self.contig.ikmers, altwindow)
        indel = VariantIndel(
            self.seqid, globalcoord - 1, refr, alt, VW=altwindow,
            RW=refrwindow, IK=str(nikmers)
        )
        yield indel

    def call_insertion(self, ksize):
        if self.targetshort:
            alt, refr, altwindow, refrwindow = self.deletion_allele(
                self.varseq, self.refrseq, ksize
            )
        else:
            refr, alt, refrwindow, altwindow = self.insertion_allele(
                self.refrseq, self.varseq, ksize
            )
        # This assertion is no longer valid when query is longer than target
        # assert len(alt) == indellength + 1
        localcoord = self.leftmatchlen
        if not self.targetshort:
            localcoord += self.offset
        globalcoord = self.cutout.local_to_global(localcoord)
        nikmers = n_ikmers_present(self.contig.ikmers, altwindow)
        indel = VariantIndel(
            self.seqid, globalcoord - 1, refr, alt, VW=altwindow,
            RW=refrwindow, IK=str(nikmers)
        )
        yield indel

    def call_indel_snvs(self, ksize, mindist=5, logstream=sys.stderr):
        offset = self.offset
        leftmatch = self.leftmatchlen
        indellength = self.indellength
        rightmatch = self.rightmatchlen

        # Left flank of the indel
        if self.targetshort:
            gdnaoffset = 0
            t = self.refrseq[:leftmatch]
            q = self.varseq[offset:offset+leftmatch]
        else:
            gdnaoffset = offset
            t = self.refrseq[offset:offset+leftmatch]
            q = self.varseq[:leftmatch]
        assert len(t) == len(q)
        diffs = [i for i in range(len(t)) if t[i] != q[i]]
        if mindist:
            diffs = trim_terminal_snvs(diffs, leftmatch, mindist, logstream)
        if len(diffs) > 0:
            for call in self.snv_variant(q, t, diffs, gdnaoffset, ksize):
                yield call

        # Right flank of the indel
        lcrf = leftmatch + indellength  # local coordinate of right flank
        if self.targetshort:
            if self.indeltype == 'I':
                gdnaoffset = leftmatch
                baseindex = offset + leftmatch + indellength
                q = self.varseq[baseindex:baseindex+rightmatch]
            else:
                gdnaoffset = leftmatch + indellength
                q = self.varseq[offset+leftmatch:offset+leftmatch+rightmatch]
            t = self.refrseq[gdnaoffset:gdnaoffset+rightmatch]
        else:
            if self.indeltype == 'I':
                gdnaoffset = offset + leftmatch
                baseindex = leftmatch + indellength
                q = self.varseq[baseindex:baseindex+rightmatch]
            else:
                gdnaoffset = offset + leftmatch + indellength
                q = self.varseq[leftmatch:leftmatch+rightmatch]
            t = self.refrseq[gdnaoffset:gdnaoffset+rightmatch]
        assert len(t) == len(q)
        diffs = [i for i in range(len(t)) if t[i] != q[i]]
        if mindist:
            diffs = trim_terminal_snvs(diffs, rightmatch, mindist, logstream)
        if len(diffs) > 0:
            for call in self.snv_variant(q, t, diffs, gdnaoffset, ksize):
                yield call


def n_ikmers_present(ikmers, window):
    n = 0
    for ikmer in ikmers:
        if ikmer.sequence in window:
            n += 1
        elif kevlar.revcom(ikmer.sequence) in window:
            n += 1
    return n


def trim_terminal_snvs(mismatches, alnlength, mindist=5, logstream=sys.stderr):
    valid = list()
    for mm in mismatches:
        if mm < mindist or alnlength - mm < mindist:
            msg = 'discarding SNV due to proximity to end of the contig'
            print('[kevlar::call] NOTE:', msg, file=logstream)
        else:
            valid.append(mm)
    return valid


def call_indel_snvs(aln, offset, ksize, leftmatch, indellength, rightmatch,
                    mindist, logstream=sys.stderr):
    calls = list()
    targetshort = False
    if offset < 0:
        offset *= -1
        targetshort = True

    # Left flank of the indel
    if targetshort:
        gdnaoffset = 0
        t = aln.refrseq[:leftmatch]
        q = aln.varseq[offset:offset+leftmatch]
    else:
        gdnaoffset = offset
        t = aln.refrseq[offset:offset+leftmatch]
        q = aln.varseq[:leftmatch]
    assert len(t) == len(q)
    diffs = [(i, t[i], q[i]) for i in range(length) if t[i] != q[i]]
    if mindist:
        diffs = trim_terminal_snvs(diffs, leftmatch, mindist, logstream)
    if len(diffs) > 0:
        snvs = snv_variant(aln, diffs, q, t, ksize)
        calls.extend(snvs)
