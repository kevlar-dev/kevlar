#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import namedtuple
import sys
import kevlar


OverlappingReadPair = namedtuple('OverlappingReadPair',
                                 'tail head offset overlap sameorient swapped')
INCOMPATIBLE_PAIR = OverlappingReadPair(None, None, None, None, None, None)


def print_read_pair(pair, position, outstream=sys.stderr):
    """Convenience print function for debugging."""
    seq2 = pair.head.sequence
    if not pair.sameorient:
        seq2 = kevlar.revcom(pair.head.sequence)
    ksize = pair.head.annotations[0].ksize

    details = '--(overlap={:d}, offset={:d}, sameorient={})-->'.format(
        pair.overlap, pair.offset, pair.sameorient
    )
    info = '[kevlar::overlap] DEBUG: shared interesting k-mer '
    info += '{:s} {:s} {:s}'.format(pair.tail.name, details, pair.head.name)

    print('≠' * 80, '\n', info, '\n', '-' * 80, '\n',
          pair.tail.sequence, '\n',
          ' ' * position, '|' * ksize, '\n',
          ' ' * pair.offset, seq2, '\n', '≠' * 80, '\n',
          sep='', file=outstream)


def check_kmer_freq_in_read_pair(read1, read2, minkmer, debugstream=None):
    """
    Check interesting k-mer frequence in each read.

    When calculating offset between a pair of reads, do not use any interesting
    k-mers that occur multiple times in either read.
    """
    maxkmer = kevlar.revcom(minkmer)
    matches1 = [k for k in read1.annotations
                if kevlar.same_seq(read1.ikmerseq(k), minkmer, maxkmer)]
    matches2 = [k for k in read2.annotations
                if kevlar.same_seq(read2.ikmerseq(k), minkmer, maxkmer)]
    nmatches1 = len(matches1)
    nmatches2 = len(matches2)
    assert nmatches1 > 0 and nmatches1 > 0, (nmatches1, nmatches2)
    if nmatches1 > 1 or nmatches2 > 1:
        if debugstream:
            message = (
                'stubbornly refusing to calculate offset bewteen {:s} and '
                '{:s}; interesting k-mer {:s} occurs multiple times'.format(
                    read1.name, read2.name, minkmer
                )
            )
            print('[kevlar::overlap] INFO', message, file=debugstream)
        return None, None

    kmer1 = matches1[0]
    kmer2 = matches2[0]
    return kmer1, kmer2


def determine_relative_orientation(read1, read2, kmer1, kmer2):
    """
    Determine the relative orientation of a pair of overlapping reads.

    Use the sequence and position of the shared interesting k-mers to determine
    the read's relative orientation.
    """
    ksize = kmer1.ksize
    pos1 = kmer1.offset
    pos2 = kmer2.offset
    kmer1seq = read1.ikmerseq(kmer1)
    kmer2seq = read2.ikmerseq(kmer2)
    sameorient = True
    if kmer1seq != kmer2seq:
        assert kmer1seq == kevlar.revcom(kmer2seq)
        sameorient = False
        if len(read1) > len(read2):
            pos1 = len(read1.sequence) - (kmer1.offset + ksize)
        else:
            pos2 = len(read2.sequence) - (kmer2.offset + ksize)

    tail, head = read1, read2
    tailpos, headpos = pos1, pos2
    read1contained = pos1 == pos2 and len(read2.sequence) > len(read1.sequence)
    if len(read1) > len(read2) or pos2 > pos1 or read1contained:
        tail, head = read2, read1
        tailpos, headpos = headpos, tailpos
    offset = tailpos - headpos

    return tail, head, offset, sameorient, tailpos


def validate_read_overlap(tail, head, offset, sameorient, minkmer, swapped):
    """Verify that the overlap between two reads is identical."""
    headseq = head.sequence if sameorient else kevlar.revcom(head.sequence)
    seg2offset = len(head.sequence) - len(tail.sequence) + offset
    if offset + len(headseq) <= len(tail.sequence):
        segment1 = tail.sequence[offset:offset+len(headseq)]
        segment2 = headseq
        seg2offset = None
    elif swapped:
        segment1 = tail.sequence[:-offset]
        segment2 = headseq[seg2offset:]
    else:
        segment1 = tail.sequence[offset:]
        segment2 = headseq[:-seg2offset]

    overlap1 = len(segment1)
    overlap2 = len(segment2)
    if overlap1 != overlap2:  # pragma: no cover
        maxkmer = kevlar.revcom(minkmer)
        print(
            '[kevlar::overlap] ERROR '
            'tail="{tail}" head="{head}" offset={offset} altoffset={altoffset}'
            ' tailoverlap={overlap} headoverlap={headover} tailolvp={tailseq}'
            ' headolvp={headseq} kmer={minkmer},{maxkmer} tailseq={tailread}'
            ' headseq={headread}'.format(
                tail=tail.name, head=tail.name,
                offset=offset,
                altoffset=seg2offset,
                overlap=overlap1, headover=len(segment2),
                tailseq=segment1, headseq=segment2,
                minkmer=minkmer, maxkmer=maxkmer,
                tailread=tail.sequence, headread=head.sequence,
            ), file=sys.stderr
        )
    assert overlap1 == overlap2
    if segment1 != segment2:
        return None
    return overlap1


def calc_offset(read1, read2, minkmer, debugstream=None):
    """
    Calculate offset between reads that share an interesting k-mer.

    Each read is annotated with its associated interesting k-mers. These are
    used to bait pairs of reads sharing interesting k-mers to build a graph of
    shared interesting k-mers. Given a pair of reads sharing an interesting,
    k-mer calculate the offset between them and determine whether they are in
    the same orientation. Any mismatches between the aligned reads (outside the
    shared k-mer) will render the offset invalid.
    """
    kmer1, kmer2 = check_kmer_freq_in_read_pair(read1, read2, minkmer,
                                                debugstream)
    if kmer1 is None or kmer2 is None:
        return INCOMPATIBLE_PAIR

    tail, head, offset, sameorient, tailpos = determine_relative_orientation(
        read1, read2, kmer1, kmer2
    )

    swapped = read1.name != tail.name and sameorient is False
    overlap = validate_read_overlap(tail, head, offset, sameorient, minkmer,
                                    swapped)
    if overlap is None:
        return INCOMPATIBLE_PAIR

    pair = OverlappingReadPair(tail=tail, head=head, offset=offset,
                               overlap=overlap, sameorient=sameorient,
                               swapped=swapped)
    if debugstream:
        print_read_pair(pair, tailpos, debugstream)

    return pair


def merge_pair(pair):
    """Assemble a pair of overlapping reads.

    Given a pair of compatible overlapping reads, collapse and merge them into
    a single sequence.
    """
    tailseq = pair.tail.sequence
    headseq = pair.head.sequence
    offset = pair.offset
    if pair.sameorient is False:
        headseq = kevlar.revcom(pair.head.sequence)
    if headseq in pair.tail.sequence:
        return pair.tail.sequence
    if pair.swapped:
        tailseq, headseq = headseq, tailseq
        offset += len(tailseq) - len(headseq)

    headindex = len(tailseq) - offset
    headsuffix = headseq[headindex:]
    tailprefix = tailseq[offset:offset+pair.overlap]
    assert tailprefix == headseq[:headindex], \
        'error: attempted to assemble incompatible reads'
    return tailseq + headsuffix
