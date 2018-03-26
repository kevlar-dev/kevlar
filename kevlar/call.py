#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from enum import Enum
import khmer
import kevlar
from kevlar.alignment import align_both_strands
import re
import sys


patterns = {
    '^(\d+)([DI])(\d+)M(\d+)[DI]$': 'snv',
    '^(\d+)([DI])(\d+)M(\d+)[DI](\d+)M$': 'snv',
    '^(\d+)([DI])(\d+)M(\d+)([ID])(\d+)M(\d+)[DI]$': 'indel',
    '^(\d+)([DI])(\d+)M(\d+)([ID])(\d+)M(\d+)[DI](\d+)M$': 'indel',
}


class VariantMapping(object):
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
        if self.alnmatch is None:
            return None
        return int(self.alnmatch.group(3))

    @property
    def indellength(self):
        if self.alnmatch is None:
            return None
        return int(self.alnmatch.group(4))

    @property
    def indeltype(self):
        if self.alnmatch is None:
            return None
        return self.alnmatch.group(5)

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
            # for call in self.call_indel_snvs(ksize):
            #     yield call
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


class Variant(object):
    """Base class for handling variant calls and no-calls."""

    def __init__(self, seqid, pos, refr, alt, **kwargs):
        """
        Constructor method.

        The `pos` parameter expects the genomic position as a 0-based index.
        Setting the `refr` or `alt` parameters to `.` will designate this
        variant as a "no call".
        """
        self._seqid = seqid
        self._pos = pos
        self._refr = refr
        self._alt = alt
        self.info = dict()
        for key, value in kwargs.items():
            self.annotate(key, value)

    @property
    def seqid(self):
        return self._seqid

    @property
    def position(self):
        return self._pos

    @property
    def vcf(self):
        """Print variant to VCF."""
        attrstr = '.'
        if len(self.info) > 0:
            kvpairs = list()
            for key in sorted(self.info):
                if key != 'QS':
                    kvpairs.append(self.attribute(key, pair=True))
            queryseq = self.attribute('QS', pair=True)
            if queryseq:
                kvpairs.append(queryseq)
            attrstr = ';'.join(kvpairs)

        filterstr = 'PASS' if self._refr != '.' else '.'
        return '{:s}\t{:d}\t.\t{:s}\t{:s}\t.\t{:s}\t{:s}'.format(
            self._seqid, self._pos + 1, self._refr, self._alt, filterstr,
            attrstr
        )

    @property
    def cigar(self):
        return self.attribute('CG')

    @property
    def window(self):
        """
        Getter method for the variant window.

        The "variant window" (abbreviated `VW` in VCF output) is the sequence
        interval in the proband contig that encompasses all k-mers overlapping
        the variant.

        GCCTAGTTAGCTAACGTCCCGATCACTGTGTCACTGC
                    .....A
                     ....A.
                      ...A..
                       ..A...
                        .A....
                         A.....
                         |        <-- position of variant
                    [---------]   <-- variant window, interval (inclusive)
                                      encompassing all 6-mers that overlap the
                                      variant
        """
        return self.attribute('VW')

    @property
    def refrwindow(self):
        """Similar to `window`, but encapsulating the reference allele."""
        return self.attribute('RW')

    def annotate(self, key, value):
        if key in self.info:
            if isinstance(self.info[key], set):
                self.info[key].add(value)
            else:
                oldvalue = self.info[key]
                self.info[key] = set((oldvalue, value))
        else:
            self.info[key] = value

    def attribute(self, key, pair=False):
        if key not in self.info:
            return None
        value = self.info[key]
        if isinstance(value, set):
            value = ','.join(sorted(value))
        value = value.replace(';', ':')
        if pair:
            keyvaluepair = '{:s}={:s}'.format(key, value)
            return keyvaluepair
        else:
            return value

    @property
    def genotypes(self):
        gt = self.attribute('GT')
        if not gt:
            return None
        return tuple(gt.split(','))


class VariantSNV(Variant):
    def __str__(self):
        return '{:s}:{:d}:{:s}->{:s}'.format(self._seqid, self._pos,
                                             self._refr, self._alt)


class VariantIndel(Variant):
    def __str__(self):
        """
        Return a string representation of this variant.

        The reason that 1 is added to the variant position is to offset the
        nucleotide shared by the reference and alternate alleles. This position
        is still 0-based (as opposed to VCF's 1-based coordinate system) but
        does not include the shared nucleotide.
        """
        pos = self._pos + 1
        if len(self._refr) > len(self._alt):
            dellength = len(self._refr) - len(self._alt)
            return '{:s}:{:d}:{:d}D'.format(self._seqid, pos, dellength)
        else:
            insertion = self._alt[1:]
            return '{:s}:{:d}:I->{:s}'.format(self._seqid, pos, insertion)


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
