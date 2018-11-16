#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from itertools import chain
import re
import sys
import kevlar
from kevlar.alignment import align_both_strands
from kevlar.cigar import AlignmentTokenizer
from kevlar.vcf import Variant
from kevlar.vcf import VariantFilter as vf


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
        self.strand = strand
        self.matedist = None
        self.trimmed = 0

        self.tok = AlignmentTokenizer(self.varseq, self.refrseq, cigar)
        self.cigar = self.tok._cigar

        self.vartype = None
        snvpattern = r'^((\d+)([DI]))?(\d+)M((\d+)[DI])?$'
        indelpattern = r'^((\d+)([DI]))?(\d+)M(\d+)([ID])(\d+)M((\d+)[DI])?$'
        if re.search(snvpattern, self.cigar):
            self.vartype = 'snv'
        elif re.search(indelpattern, self.cigar):
            self.vartype = 'indel'

    def __str__(self):
        fulltarget, fullquery = '', ''
        for token in self.tok.blocks:
            fulltarget += token.target if token.target else '-' * token.length
            fullquery += token.query if token.query else '-' * token.length

        fullmatch = list()
        for t, q in zip(fulltarget, fullquery):
            c = '|' if t == q else ' '
            fullmatch.append(c)
        fullmatch = ''.join(fullmatch)

        outlines = list()
        i = 0
        while i < len(fulltarget):
            outlines.append(fulltarget[i:i+80])
            outlines.append(fullmatch[i:i+80])
            outlines.append(fullquery[i:i+80])
            outlines.append('')
            i += 80
        return '\n'.join(outlines).strip()

    @property
    def interval(self):
        return self.cutout.interval

    @property
    def ikmers(self):
        for kmer in self.contig.annotations:
            seq = self.contig.ikmerseq(kmer)
            yield seq
            yield kevlar.revcom(seq)

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
        if self.vartype is None:
            return None
        if self.tok.blocks[0].type == 'M':
            return 0
        return self.tok.blocks[0].length

    @property
    def targetshort(self):
        if self.vartype is None:
            return None
        return self.tok.blocks[0].type == 'I'

    @property
    def match(self):
        if self.vartype != 'snv':
            return None
        i = 0 if self.tok.blocks[0].type == 'M' else 1
        return self.tok.blocks[i]

    @property
    def leftflank(self):
        if self.vartype != 'indel':
            return None
        i = 0 if self.tok.blocks[0].type == 'M' else 1
        return self.tok.blocks[i]

    @property
    def indel(self):
        if self.vartype != 'indel':
            return None
        i = 1 if self.tok.blocks[0].type == 'M' else 2
        return self.tok.blocks[i]

    @property
    def indeltype(self):
        if self.vartype != 'indel':
            return None
        return self.indel.type

    @property
    def rightflank(self):
        if self.vartype != 'indel':
            return None
        i = -1 if self.tok.blocks[-1].type == 'M' else -2
        return self.tok.blocks[i]

    def is_passenger(self, call):
        if call.window is None:
            return False
        numikmers = sum([1 for k in self.ikmers if k in call.window])
        return numikmers == 0

    def call_variants(self, ksize, mindist=6, logstream=sys.stderr):
        """Attempt to call variants from this contig alignment.

        If the alignment CIGAR matches a known pattern, the appropriate caller
        is invoked (SNV or INDEL caller). If not, a "no call" is reported.

        If an SNV call is within `mindist` base pairs of the end of the
        alignment it is ignored. Set to `None` to disable this behavior.

        Variant calls with no spanning interesting k-mers are designated as
        "passenger calls" and discarded.
        """
        offset = 0 if self.targetshort else self.offset
        if self.vartype == 'snv':
            caller = self.call_snv(self.match.query, self.match.target, offset,
                                   ksize, mindist, logstream=logstream)
            for call in caller:
                if self.is_passenger(call):
                    call.filter(vf.PassengerVariant)
                yield call
        elif self.vartype == 'indel':
            indelcaller = self.call_indel(ksize)
            indel = next(indelcaller)
            if self.is_passenger(indel):
                indel.filter(vf.PassengerVariant)
            yield indel

            leftflankcaller = self.call_snv(
                self.leftflank.query, self.leftflank.target, offset, ksize,
                mindist, donocall=False
            )
            offset += self.leftflank.length
            if self.indeltype == 'D':
                offset += self.indel.length
            rightflankcaller = self.call_snv(
                self.rightflank.query, self.rightflank.target, offset, ksize,
                mindist, donocall=False
            )
            for call in chain(leftflankcaller, rightflankcaller):
                if self.is_passenger(call):
                    call.filter(vf.PassengerVariant)
                yield call
        else:
            nocall = Variant(
                self.seqid, self.pos, '.', '.', CONTIG=self.varseq,
                CIGAR=self.cigar, KSW2=str(self.score)
            )
            nocall.filter(vf.InscrutableCigar)
            yield nocall

    def call_snv(self, qseq, tseq, offset, ksize, mindist=6, donocall=True,
                 logstream=sys.stderr):
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
        if length < ksize:
            return
        diffs = [i for i in range(length) if tseq[i] != qseq[i]]
        if mindist:
            self.trimmed, diffs = trim_terminal_snvs(diffs, length, mindist)
        if len(diffs) == 0 or len(diffs) > 4:
            if donocall:
                nocall = Variant(
                    self.seqid, self.cutout.local_to_global(offset), '.', '.',
                    CONTIG=qseq, CIGAR=self.cigar, KSW2=str(self.score),
                    IKMERS=str(len(self.contig.annotations))
                )
                if len(diffs) == 0:
                    nocall.filter(vf.PerfectMatch)
                if len(diffs) > 4:
                    nocall.filter(vf.NumerousMismatches)
                yield nocall
            return

        for pos in diffs:
            minpos = max(pos - ksize + 1, 0)
            maxpos = min(pos + ksize, length)
            altwindow = qseq[minpos:maxpos]
            refrwindow = tseq[minpos:maxpos]

            refr = tseq[pos].upper()
            alt = qseq[pos].upper()
            localcoord = pos + offset
            globalcoord = self.cutout.local_to_global(localcoord)
            nikmers = n_ikmers_present(self.contig, altwindow)
            snv = Variant(
                self.seqid, globalcoord, refr, alt, CONTIG=qseq,
                CIGAR=self.cigar, KSW2=str(self.score), IKMERS=str(nikmers),
                ALTWINDOW=altwindow, REFRWINDOW=refrwindow
            )
            yield snv

    def call_indel(self, ksize):
        if self.indeltype == 'D':
            refrwindow = self.leftflank.target[-(ksize-1):] \
                + self.indel.target \
                + self.rightflank.target[:(ksize-1)]
            refrallele = self.leftflank.target[-1] + self.indel.target
            altwindow = self.leftflank.query[-(ksize-1):] \
                + self.rightflank.query[:(ksize-1)]
            altallele = self.leftflank.query[-1]
        else:
            refrwindow = self.leftflank.target[-(ksize-1):] \
                + self.rightflank.target[:(ksize-1)]
            refrallele = self.leftflank.target[-1]
            altwindow = self.leftflank.query[-(ksize-1):] \
                + self.indel.query \
                + self.rightflank.query[:(ksize-1)]
            altallele = self.leftflank.query[-1] + self.indel.query
        nikmers = n_ikmers_present(self.contig, altwindow)
        localcoord = 0 if self.targetshort else self.offset
        localcoord += self.leftflank.length
        globalcoord = self.cutout.local_to_global(localcoord)
        indel = Variant(
            self.seqid, globalcoord - 1, refrallele, altallele,
            CONTIG=self.refrseq, CIGAR=self.cigar, KSW2=str(self.score),
            IKMERS=str(nikmers), ALTWINDOW=altwindow, REFRWINDOW=refrwindow
        )
        yield indel


def n_ikmers_present(record, window):
    n = 0
    for ikmer in record.annotations:
        seq = record.ikmerseq(ikmer)
        if seq in window:
            n += 1
        elif kevlar.revcom(seq) in window:
            n += 1
    return n


def trim_terminal_snvs(mismatches, alnlength, mindist=5, logstream=sys.stderr):
    valid = list()
    trimcount = 0
    for mm in mismatches:
        if mm < mindist or alnlength - mm < mindist:
            trimcount += 1
            # msg = 'discarding SNV due to proximity to end of the contig'
            # print('[kevlar::call] NOTE:', msg, file=logstream)
        else:
            valid.append(mm)
    return trimcount, valid
