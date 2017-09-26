#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import re
import sys
import khmer
import kevlar


class Variant(object):
    def __init__(self, seqid, pos, refr, alt, kmers=None):
        self._seqid = seqid
        self._pos = pos
        self._refr = refr
        self._alt = alt
        self.kmers = kmers
        self.info = dict()
        if kmers:
            kmerstr = ','.join(kmers)
            self.info['KevlarKmers'] = kmerstr

    @property
    def vcf(self):
        info = '.'
        if len(self.info) > 0:
            infokvp = [
                '{:s}={:s}'.format(key, value.replace(';', ':'))
                for key, value in self.info.items()
            ]
            info = ';'.join(infokvp)

        return '{:s}\t{:d}\t.\t{:s}\t{:s}\t.\tPASS\t{:s}'.format(
            self._seqid, self._pos + 1, self._refr, self._alt, info
        )


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


def local_to_global(localcoord, subseqid):
    match = re.search('(\S+)_(\d+)-(\d+)', subseqid)
    assert match, 'unable to parse subseqid {:s}'.format(subseqid)
    seqid = match.group(1)
    globaloffset = int(match.group(2))
    globalcoord = globaloffset + localcoord
    return seqid, globalcoord


def call_snv(target, query, offset, length, ksize):
    t = target.sequence[offset:offset+length]
    q = query.sequence[:length]
    diffs = [(i, t[i], q[i]) for i in range(length) if t[i] != q[i]]
    if len(diffs) == 0:
        raise StopIteration

    snvs = list()
    for diff in diffs:
        minpos = max(diff[0] - ksize, 0)
        maxpos = min(diff[0] + ksize, length)
        window = t[minpos:maxpos]
        kmers = [t[i:i+ksize] for i in range(ksize)]

        refr = diff[1].upper()
        alt = diff[2].upper()
        localcoord = offset + diff[0]
        seqid, globalcoord = local_to_global(localcoord, target.name)
        snv = VariantSNV(seqid, globalcoord, refr, alt, kmers=kmers)
        snvs.append(snv)
    return snvs


def call_deletion(target, query, offset, leftmatch, indellength):
    refr = target.sequence[offset+leftmatch-1:offset+leftmatch+indellength]
    alt = refr[0]
    assert len(refr) == indellength + 1
    localcoord = offset + leftmatch
    seqid, globalcoord = local_to_global(localcoord, target.name)
    return [VariantIndel(seqid, globalcoord - 1, refr, alt)]


def call_insertion(target, query, offset, leftmatch, indellength):
    insertion = query.sequence[leftmatch-1:leftmatch+indellength]
    refr = insertion[0]
    assert len(insertion) == indellength + 1
    localcoord = offset + leftmatch
    seqid, globalcoord = local_to_global(localcoord, target.name)
    return [VariantIndel(seqid, globalcoord - 1, refr, insertion)]


def make_call(target, query, cigar, ksize):
    snvmatch = re.search('^(\d+)D(\d+)M(\d+)D$', cigar)
    snvmatch2 = re.search('^(\d+)D(\d+)M(\d+)D(\d+)M$', cigar)
    if snvmatch:
        offset = int(snvmatch.group(1))
        length = int(snvmatch.group(2))
        return call_snv(target, query, offset, length, ksize)
    elif snvmatch2 and int(snvmatch2.group(4)) <= 5:
        offset = int(snvmatch2.group(1))
        length = int(snvmatch2.group(2))
        return call_snv(target, query, offset, length, ksize)

    indelmatch = re.search('^(\d+)D(\d+)M(\d+)([ID])(\d+)M(\d+)D$', cigar)
    indelmatch2 = re.search('^(\d+)D(\d+)M(\d+)([ID])(\d+)M(\d+)D(\d+)M$',
                            cigar)
    if indelmatch:
        offset = int(indelmatch.group(1))
        leftmatch = int(indelmatch.group(2))
        indellength = int(indelmatch.group(3))
        indeltype = indelmatch.group(4)
        callfunc = call_deletion if indeltype == 'D' else call_insertion
        return callfunc(target, query, offset, leftmatch, indellength)
    elif indelmatch2 and int(indelmatch2.group(7)) <= 5:
        offset = int(indelmatch2.group(1))
        leftmatch = int(indelmatch2.group(2))
        indellength = int(indelmatch2.group(3))
        indeltype = indelmatch2.group(4)
        callfunc = call_deletion if indeltype == 'D' else call_insertion
        return callfunc(target, query, offset, leftmatch, indellength)


def call(targetlist, querylist, match=1, mismatch=2, gapopen=5, gapextend=0,
         ksize=31):
    """
    Wrap the `kevlar call` procedure as a generator function.

    Input is the following.
    - an iterable containing one or more target sequences from the reference
      genome, stored as khmer or screed sequence records
    - an iterable containing one or more contigs assembled by kevlar, stored as
      khmer or screed sequence records
    - alignment match score (integer)
    - alignment mismatch penalty (integer)
    - alignment gap open penalty (integer)
    - alignment gap extension penalty (integer)

    The function yields tuples of target sequence name, query sequence name,
    and alignment CIGAR string
    """
    for target in sorted(targetlist, key=lambda record: record.name):
        for query in sorted(querylist, reverse=True, key=len):
            cigar = kevlar.align(target.sequence, query.sequence, match,
                                 mismatch, gapopen, gapextend)
            varcalls = make_call(target, query, cigar, ksize)
            if varcalls is not None:
                for varcall in varcalls:
                    varcall.info['query'] = query.name
                    varcall.info['CIGAR'] = cigar
                    yield target.name, query.name, cigar, varcall
            else:
                yield target.name, query.name, cigar, None


def main(args):
    outstream = kevlar.open(args.out, 'w')
    qinstream = kevlar.parse_augmented_fastx(kevlar.open(args.queryseq, 'r'))
    queryseqs = [record for record in qinstream]
    targetseqs = [record for record in khmer.ReadParser(args.targetseq)]
    caller = call(
        targetseqs, queryseqs,
        args.match, args.mismatch, args.open, args.extend
    )
    for target, query, cigar, varcall in caller:
        if varcall is None:
            print('#', target, query, cigar, file=outstream)
        else:
            print(varcall.vcf, file=outstream)
