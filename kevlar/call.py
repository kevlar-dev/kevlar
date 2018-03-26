#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import khmer
import kevlar
import re
import sys


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

    def call_variants(self, ksize, mindist=5, logstream=sys.stderr):
        snvmatch = re.search('^(\d+)([DI])(\d+)M(\d+)[DI]$', self.cigar)
        snvmatch2 = re.search('^(\d+)([DI])(\d+)M(\d+)[DI](\d+)M$', self.cigar)
        if snvmatch:
            offset = int(snvmatch.group(1))
            if snvmatch.group(2) == 'I':
                offset *= -1
            length = int(snvmatch.group(3))
            return call_snv(self, offset, length, ksize, mindist, logstream)
        elif snvmatch2 and int(snvmatch2.group(5)) <= 5:
            offset = int(snvmatch2.group(1))
            if snvmatch2.group(2) == 'I':
                offset *= -1
            length = int(snvmatch2.group(3))
            return call_snv(self, offset, length, ksize, mindist, logstream)

        indelmatch = re.search(
            '^(\d+)([DI])(\d+)M(\d+)([ID])(\d+)M(\d+)[DI]$', self.cigar
        )
        indelmatch2 = re.search(
            '^(\d+)([DI])(\d+)M(\d+)([ID])(\d+)M(\d+)[DI](\d+)M$', self.cigar
        )
        indmatch = None
        if indelmatch:
            indmatch = indelmatch
        elif indelmatch2 and int(indelmatch2.group(8)) <= 5:
            indmatch = indelmatch2
        if indmatch:
            offset = int(indmatch.group(1))
            if indmatch.group(2) == 'I':
                offset *= -1
            leftmatchlen = int(indmatch.group(3))
            indellength = int(indmatch.group(4))
            indeltype = indmatch.group(5)
            rightmatchlen = int(indmatch.group(6))
            callfunc = call_deletion if indeltype == 'D' else call_insertion
            indels = callfunc(self, offset, ksize, leftmatchlen, indellength)
            snvs = call_indel_snvs(self, offset, ksize, leftmatchlen,
                                   rightmatchlen)
            return indels + snvs

        nocall = Variant(
            self.seqid, self.pos, '.', '.', NC='inscrutablecigar',
            QN=self.contig.name, QS=self.contig.sequence, CG=self.cigar,
        )
        return [nocall]


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


def snv_variant(aln, mismatches, query, target, ksize):
    snvs = list()
    length = len(query)
    for diff in mismatches:
        minpos = max(diff[0] - ksize + 1, 0)
        maxpos = min(diff[0] + ksize, length)
        window = query[minpos:maxpos]
        refrwindow = target[minpos:maxpos]

        # numoverlappingkmers = len(window) - ksize + 1
        # kmers = [window[i:i+ksize] for i in range(numoverlappingkmers)]
        refr = diff[1].upper()
        alt = diff[2].upper()
        localcoord = diff[0]
        if not targetshort:
            localcoord += offset
        globalcoord = aln.cutout.local_to_global(localcoord)
        nikmers = n_ikmers_present(aln.contig.ikmers, window)
        snv = VariantSNV(aln.seqid, globalcoord, refr, alt, VW=window,
                         RW=refrwindow, IK=str(nikmers))
        snvs.append(snv)
    return snvs


def trim_terminal_snvs(mismatches, alnlength, mindist=5, logstream=sys.stderr):
    valid = list()
    for mm in mismatches:
        if mm[0] < mindist or alnlength - mm[0] < mindist:
            msg = 'discarding SNV due to proximity to end of the contig'
            print('[kevlar::call] NOTE:', msg, file=logstream)
        else:
            valid.append(d)
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


def call_snv(aln, offset, length, ksize, mindist, logstream=sys.stderr):
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
    if mindist:
        diffs = trim_terminal_snvs(diffs, length, mindist, logstream)
    if len(diffs) == 0:
        nocall = Variant(aln.seqid, aln.cutout.local_to_global(gdnaoffset),
                         '.', '.', NC='perfectmatch', QN=aln.contig.name, QS=q)
        return [nocall]

    snvs = snv_variant(aln, diffs, q, t, ksize)
    return snvs


def deletion_allele(target, query, offset, ksize, leftmatch, indellength):
    minpos = max(leftmatch - ksize + 1, 0)
    maxpos = min(leftmatch + ksize - 1, len(query))
    altwindow = query[minpos:maxpos]
    minpos += offset
    maxpos += offset + indellength
    refrwindow = target[minpos:maxpos]

    refr = target[offset+leftmatch-1:offset+leftmatch+indellength]
    alt = refr[0]
    return refr, alt, refrwindow, altwindow


def insertion_allele(target, query, offset, ksize, leftmatch, indellength):
    minpos = max(leftmatch - ksize + 1, 0)
    maxpos = min(leftmatch + ksize + indellength - 1, len(query))
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
    nikmers = n_ikmers_present(aln.contig.ikmers, altwindow)
    return [VariantIndel(aln.seqid, globalcoord - 1, refr, alt, VW=altwindow,
                         RW=refrwindow, IK=str(nikmers))]


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
    nikmers = n_ikmers_present(aln.contig.ikmers, altwindow)
    return [VariantIndel(aln.seqid, globalcoord - 1, refr, alt, VW=altwindow,
                         RW=refrwindow, IK=str(nikmers))]


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
