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


def local_to_global(localcoord, subseqid):
    match = re.search('(\S+)_(\d+)-(\d+)', subseqid)
    assert match, 'unable to parse subseqid {:s}'.format(subseqid)
    seqid = match.group(1)
    globaloffset = int(match.group(2))
    globalcoord = globaloffset + localcoord
    return seqid, globalcoord


def call_snv(target, query, offset, length, ksize):
    targetshort = False
    if offset < 0:
        offset *= -1
        gdnaoffset = 0
        targetshort = True
        t = target.sequence[:length]
        q = query.sequence[offset:offset+length]
    else:
        gdnaoffset = offset
        t = target.sequence[offset:offset+length]
        q = query.sequence[:length]
    diffs = [(i, t[i], q[i]) for i in range(length) if t[i] != q[i]]
    if len(diffs) == 0:
        seqid, globalcoord = local_to_global(gdnaoffset, target.name)
        nocall = Variant(seqid, globalcoord, '.', '.', NC='perfectmatch',
                         QN=query.name, QS=q)
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
        seqid, globalcoord = local_to_global(localcoord, target.name)
        snv = Variant(seqid, globalcoord, refr, alt, VW=window,
                      RW=refrwindow, IK=str(len(query.ikmers)))
        snvs.append(snv)
    return snvs


def deletion_allele(target, query, offset, ksize, leftmatch, indellength):
    minpos = leftmatch - ksize + 1
    maxpos = leftmatch + ksize - 1
    altwindow = query.sequence[minpos:maxpos]
    minpos += offset
    maxpos += offset + indellength
    refrwindow = target.sequence[minpos:maxpos]

    refr = target.sequence[offset+leftmatch-1:offset+leftmatch+indellength]
    alt = refr[0]
    return refr, alt, refrwindow, altwindow


def insertion_allele(target, query, offset, ksize, leftmatch, indellength):
    minpos = leftmatch - ksize + 1
    maxpos = leftmatch + ksize + indellength - 1
    altwindow = query.sequence[minpos:maxpos]
    minpos += offset
    maxpos += offset - indellength
    refrwindow = target.sequence[minpos:maxpos]

    alt = query.sequence[leftmatch-1:leftmatch+indellength]
    refr = alt[0]
    return refr, alt, refrwindow, altwindow


def call_deletion(target, query, offset, ksize, leftmatch, indellength):
    if offset < 0:
        offset *= -1
        targetshort = True
        alt, refr, altwindow, refrwindow = insertion_allele(
            query, target, offset, ksize, leftmatch, indellength
        )
    else:
        targetshort = False
        refr, alt, refrwindow, altwindow = deletion_allele(
            target, query, offset, ksize, leftmatch, indellength
        )
    # This assertion is no longer valid when query is longer than target
    # assert len(refr) == indellength + 1
    localcoord = leftmatch
    if not targetshort:
        localcoord += offset
    seqid, globalcoord = local_to_global(localcoord, target.name)
    var = Variant(seqid, globalcoord - 1, refr, alt, VW=altwindow,
                  RW=refrwindow, IK=str(len(query.ikmers)))
    return [var]


def call_insertion(target, query, offset, ksize, leftmatch, indellength):
    if offset < 0:
        offset *= -1
        targetshort = True
        alt, refr, altwindow, refrwindow = deletion_allele(
            query, target, offset, ksize, leftmatch, indellength
        )
    else:
        targetshort = False
        refr, alt, refrwindow, altwindow = insertion_allele(
            target, query, offset, ksize, leftmatch, indellength
        )

    assert len(alt) == indellength + 1
    localcoord = leftmatch
    if not targetshort:
        localcoord += offset
    seqid, globalcoord = local_to_global(localcoord, target.name)
    var = Variant(seqid, globalcoord - 1, refr, alt, VW=altwindow,
                  RW=refrwindow, IK=str(len(query.ikmers)))
    return [var]


def make_call(target, query, cigar, ksize):
    snvmatch = re.search('^(\d+)([DI])(\d+)M(\d+)[DI]$', cigar)
    snvmatch2 = re.search('^(\d+)([DI])(\d+)M(\d+)[DI](\d+)M$', cigar)
    if snvmatch:
        offset = int(snvmatch.group(1))
        if snvmatch.group(2) == 'I':
            offset *= -1
        length = int(snvmatch.group(3))
        return call_snv(target, query, offset, length, ksize)
    elif snvmatch2 and int(snvmatch2.group(5)) <= 5:
        offset = int(snvmatch2.group(1))
        if snvmatch2.group(2) == 'I':
            offset *= -1
        length = int(snvmatch2.group(3))
        return call_snv(target, query, offset, length, ksize)

    indelmatch = re.search(
        '^(\d+)([DI])(\d+)M(\d+)([ID])(\d+)M(\d+)[DI]$', cigar
    )
    indelmatch2 = re.search(
        '^(\d+)([DI])(\d+)M(\d+)([ID])(\d+)M(\d+)[DI](\d+)M$', cigar
    )
    if indelmatch:
        offset = int(indelmatch.group(1))
        if indelmatch.group(2) == 'I':
            offset *= -1
        leftmatch = int(indelmatch.group(3))
        indellength = int(indelmatch.group(4))
        indeltype = indelmatch.group(5)
        callfunc = call_deletion if indeltype == 'D' else call_insertion
        return callfunc(target, query, offset, ksize, leftmatch, indellength)
    elif indelmatch2 and int(indelmatch2.group(8)) <= 5:
        offset = int(indelmatch2.group(1))
        if indelmatch2.group(2) == 'I':
            offset *= -1
        leftmatch = int(indelmatch2.group(3))
        indellength = int(indelmatch2.group(4))
        indeltype = indelmatch2.group(5)
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


def call(targetlist, querylist, match=1, mismatch=2, gapopen=5, gapextend=0,
         ksize=31, casecounts=None, controlcounts=None, refr=None, mu=30.0,
         sigma=8.0, epsilon=0.01, caselabel=None, ctrllabels=None):
    """
    Wrap the `kevlar call` procedure as a generator function.

    The `targetlist` and `querylist` variables should be iterables containig
    sequences to be aligned. These lists will be sorted and iterated over
    multiple times, so they cannot be used in a streaming fashion.

    The `match, `mismatch`, `gapopen`, and `gapextend` parameters are for
    alignment scoring and should all be non-negative integers.

    If `casecounts` and `controlcounts` are provided, the likelihood of each
    variant will be computed. `casecounts` should be a Counttable with k-mer
    abundances for the proband sample, and `controlcounts` should be a list of
    Counttables, one per parent/control sample. If `refr` is provided (a sketch
    with k-mers from the reference genome), any k-mers containing the alternate
    allele that occur elsewhere in the genome will be discarded.

    The `mu` and `sigma` parameters refer to the mean and standard deviation of
    the observed k-mer abundance distributions, and each can be specified using
    a single float or a list of floats (1 per sample). The `epsilon` parameter
    is an error rate, and can also be specified using a float or a list of
    floats.
    """
    dolike = casecounts is not None and controlcounts is not None
    varcalls = list()
    for query in sorted(querylist, reverse=True, key=len):
        alignments = list()
        for target in sorted(targetlist, key=lambda record: record.name):
            cigar, score, strand = align_both_strands(
                target.sequence, query.sequence, match, mismatch, gapopen,
                gapextend
            )
            alignments.append((target, cigar, score, strand))
        alignments.sort(key=lambda a: a[2], reverse=True)
        if len(alignments) == 1:
            aligns2report = alignments
        else:
            scrtbl = [a for a in alignments if alignment_interpretable(a[1])]
            if len(scrtbl) == 0:
                finallist = alignments
            else:
                finallist = scrtbl
            bestscore = finallist[0][2]
            aligns2report = [a for a in finallist if a[2] == bestscore]

        for alignment in aligns2report:
            besttarget, bestcigar, bestscore, bestorientation = alignment
            if bestorientation == -1:
                query.sequence = kevlar.revcom(query.sequence)
            for varcall in make_call(besttarget, query, bestcigar, ksize):
                # if dolike:
                #     compute_likelihoods(
                #         variant, casecounts, controlcounts, refr=refr,
                #         mean=mu, sd=sigma, error=epsilon,
                #         caselabel=caselabel, ctrllables=ctrllabels
                #     )
                # varcalls.append(varcall)
                yield varcall
            if bestorientation == -1:
                # Change it back!
                # There's a better way to do this, but this works for now.
                query.sequence = kevlar.revcom(query.sequence)

    # if dolike:
    #     varcalls.sort(key=lambda v: float(v.info['DN']), reverse=True)
    # for v in varcalls:
    #     yield v


def main(args):
    outstream = kevlar.open(args.out, 'w')
    writer = kevlar.vcf.VCFWriter(source='kevlar::call')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(args.queryseq, 'r'))
    queryseqs = list(qinstream)
    targetseqs = list(khmer.ReadParser(args.targetseq))

    caselabel = None
    ctrllabels = None
    if args.like_filter:
        refr = khmer.Nodetable.load(args.refr)
        casecounts = khmer.Counttable.load(args.case)
        ctrlcounts = [khmer.Counttable.load(c) for c in args.control]
        caselabel = args.case_label if args.case_label else 'Case'
        writer.register_sample(caselabel)
        if args.ctrl_labels:
            ctrllabels = args.ctrl_labels.split(',')
        else:
            numcontrols = len(args.control)
            ctrllabels = ['Control' + str(i+1) for i in range(numcontrols)]
        for label in ctrllabels:
            writer.register_sample(label)

    caller = call(
        targetseqs, queryseqs, args.match, args.mismatch, args.open,
        args.extend, args.ksize, args.case, args.ctrl, args.refr, args.mu,
        args.sigma, args.epsilon, caselabel, ctrllabels
    )

    for varcall in caller:
        writer.write(varcall, fh=outstream)
