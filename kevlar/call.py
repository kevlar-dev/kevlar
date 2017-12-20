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
from kevlar.vcf import Variant


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
        seqid, globalcoord = local_to_global(offset, target.name)
        nocall = Variant(seqid, globalcoord, '.', '.', NC='perfectmatch', CS=q)
        return [nocall]

    snvs = list()
    for diff in diffs:
        minpos = max(diff[0] - ksize + 1, 0)
        maxpos = min(diff[0] + ksize, length)
        window = q[minpos:maxpos]
        refrwindow = t[minpos:maxpos]

        numoverlappingkmers = len(window) - ksize + 1
        kmers = [window[i:i+ksize] for i in range(numoverlappingkmers)]
        refr = diff[1].upper()
        alt = diff[2].upper()
        localcoord = offset + diff[0]
        seqid, globalcoord = local_to_global(localcoord, target.name)
        snv = Variant(seqid, globalcoord, refr, alt, VW=window, RW=refrwindow)
        snvs.append(snv)
    return snvs


def call_deletion(target, query, offset, ksize, leftmatch, indellength):
    minpos = leftmatch - ksize + 1
    maxpos = leftmatch + ksize - 1
    window = query.sequence[minpos:maxpos]
    minpos += offset
    maxpos += offset + indellength
    refrwindow = target.sequence[minpos:maxpos]

    refr = target.sequence[offset+leftmatch-1:offset+leftmatch+indellength]
    alt = refr[0]
    assert len(refr) == indellength + 1
    localcoord = offset + leftmatch
    seqid, globalcoord = local_to_global(localcoord, target.name)
    var = Variant(seqid, globalcoord - 1, refr, alt, VW=window, RW=refrwindow)
    return [var]


def call_insertion(target, query, offset, ksize, leftmatch, indellength):
    minpos = leftmatch - ksize + 1
    maxpos = leftmatch + ksize + indellength - 1
    window = query.sequence[minpos:maxpos]
    minpos += offset
    maxpos += offset - indellength
    refrwindow = target.sequence[minpos:maxpos]

    alt = query.sequence[leftmatch-1:leftmatch+indellength]
    refr = alt[0]
    insertion = alt[1:]
    assert len(insertion) == indellength
    localcoord = offset + leftmatch
    seqid, globalcoord = local_to_global(localcoord, target.name)
    var = Variant(seqid, globalcoord - 1, refr, alt, VW=window, RW=refrwindow)
    return [var]


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
        return callfunc(target, query, offset, ksize, leftmatch, indellength)
    elif indelmatch2 and int(indelmatch2.group(7)) <= 5:
        offset = int(indelmatch2.group(1))
        leftmatch = int(indelmatch2.group(2))
        indellength = int(indelmatch2.group(3))
        indeltype = indelmatch2.group(4)
        callfunc = call_deletion if indeltype == 'D' else call_insertion
        return callfunc(target, query, offset, ksize, leftmatch, indellength)

    seqid, globalcoord = local_to_global(0, target.name)
    nocall = Variant(seqid, globalcoord, '.', '.', NC='inscrutablecigar',
                     CS=query.sequence, CIGAR=cigar)
    return [nocall]


def align_both_strands(targetseq, queryseq, match=1, mismatch=2, gapopen=5,
                       gapextend=0):
    cigar1, score1 = kevlar.align(targetseq, queryseq, match, mismatch,
                                  gapopen, gapextend)
    cigar2, score2 = kevlar.align(targetseq, kevlar.revcom(queryseq), match,
                                  mismatch, gapopen, gapextend)

    if score2 > score1:
        cigar = cigar2
        score = score2
        strand = -1
    else:
        cigar = cigar1
        score = score1
        strand = 1
    return cigar, score, strand


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
    for query in sorted(querylist, reverse=True, key=len):
        bestcigar = None
        bestscore = None
        besttarget = None
        bestorientation = None
        for target in sorted(targetlist, key=lambda record: record.name):
            cigar, score, strand = align_both_strands(
                target.sequence, query.sequence, match, mismatch, gapopen,
                gapextend
            )
            if bestscore is None or score > bestscore:
                bestscore = score
                bestcigar = cigar
                besttarget = target
                bestorientation = strand

        if bestorientation == -1:
            query.sequence = kevlar.revcom(query.sequence)
        for varcall in make_call(besttarget, query, bestcigar, ksize):
            yield varcall


def main(args):
    outstream = kevlar.open(args.out, 'w')
    qinstream = kevlar.parse_augmented_fastx(kevlar.open(args.queryseq, 'r'))
    queryseqs = [record for record in qinstream]
    targetseqs = [record for record in khmer.ReadParser(args.targetseq)]
    caller = call(
        targetseqs, queryseqs,
        args.match, args.mismatch, args.open, args.extend,
        args.ksize,
    )
    for varcall in caller:
        print(varcall.vcf, file=outstream)
