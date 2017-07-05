#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict, namedtuple
import itertools
import sys
import networkx
import screed
import kevlar


OverlappingReadPair = namedtuple('OverlappingReadPair',
                                 'tail head offset overlap sameorient swapped')
INCOMPATIBLE_PAIR = OverlappingReadPair(None, None, None, None, None, None)


def print_read_pair(pair, position, outstream=sys.stderr):
    """Convenience print function for debugging."""
    seq2 = pair.head.sequence
    if not pair.sameorient:
        seq2 = kevlar.revcom(pair.head.sequence)
    ksize = len(pair.head.ikmers[0].sequence)

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
    matches1 = [k for k in read1.ikmers
                if kevlar.same_seq(k.sequence, minkmer, maxkmer)]
    matches2 = [k for k in read2.ikmers
                if kevlar.same_seq(k.sequence, minkmer, maxkmer)]
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
    ksize = len(kmer1.sequence)
    pos1 = kmer1.offset
    pos2 = kmer2.offset
    sameorient = True
    if kmer1.sequence != kmer2.sequence:
        assert kmer1.sequence == kevlar.revcom(kmer2.sequence)
        sameorient = False
        pos2 = len(read2.sequence) - (kmer2.offset + ksize)

    tail, head = read1, read2
    tailpos, headpos = pos1, pos2
    read1contained = pos1 == pos2 and len(read2.sequence) > len(read1.sequence)
    if pos2 > pos1 or read1contained:
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


def graph_add_edge(graph, pair, minkmer):
    """
    Add edge between two nodes in the "shared interesting k-mer" read graph.

    If the edge already exists, make sure that the existing edge matches the
    edge that would have been added.
    """
    tailname, headname = pair.tail.name, pair.head.name
    if tailname in graph and headname in graph[tailname]:
        assert graph[tailname][headname]['offset'] == pair.offset
        if graph[tailname][headname]['tail'] == tailname:
            assert graph[tailname][headname]['overlap'] == pair.overlap
        graph[tailname][headname]['ikmers'].add(minkmer)
    else:
        graph.add_edge(tailname, headname, offset=pair.offset,
                       overlap=pair.overlap, ikmers=set([minkmer]),
                       orient=pair.sameorient, tail=tailname,
                       swapped=pair.swapped)


def graph_init_abund_check_pass(numreads, minkmer, minabund=5, maxabund=500,
                                logstream=None):
    """Check whether the k-mer falls within the expected range of abundance."""
    if maxabund and numreads > maxabund:
        msg = '            skipping k-mer with abundance {:d}'.format(numreads)
        print(msg, file=logstream)
        return False
    if minabund and numreads < minabund:
        message = '[kevlar::overlap] WARNING: k-mer {}'.format(minkmer)
        message += ' (rev. comp. {})'.format(kevlar.revcom(minkmer))
        message += ' only has abundance {:d}'.format(numreads)
        out = logstream if logstream is not None else sys.stderr
        print(message, file=out)
        return False
    return True


def graph_init_strict(reads, kmers, minabund=5, maxabund=500, logstream=None):
    """
    Initialize the "shared interesting k-mer" read graph.

    Iterate through each interesting k-mer, consider every pair of reads
    containing that k-mer, and determine whether there should be a edge between
    the pair of reads in the graph.
    """
    graph = networkx.Graph()
    nkmers = len(kmers)
    for n, minkmer in enumerate(kmers, 1):
        if n % 100 == 0 and logstream:  # pragma: no cover
            msg = 'processed {:d}/{:d} shared novel k-mers'.format(n, nkmers)
            print('[kevlar::overlap]    ', msg, sep='', file=logstream)

        readnames = kmers[minkmer]
        if not graph_init_abund_check_pass(len(readnames), minkmer, minabund,
                                           maxabund, logstream):
            continue

        readset = [reads[rn] for rn in readnames]
        for read1, read2 in itertools.combinations(readset, 2):
            pair = kevlar.overlap.calc_offset(read1, read2, minkmer, logstream)
            if pair is kevlar.overlap.INCOMPATIBLE_PAIR:
                # Shared k-mer but bad overlap
                continue
            graph_add_edge(graph, pair, minkmer)
    return graph


def graph_init_basic(reads_by_kmer, upint=1000, logstream=None):
    """
    Initialize read graph by shared interesting k-mers.

    Reads do not require perfect matching in the overlap to be connected by an
    edge in the graph.
    """
    read_graph = networkx.Graph()
    for n, kmer in enumerate(reads_by_kmer):
        if logstream and n > 0 and n % upint == 0:
            print('    build shared novel k-mer graph:', n, file=logstream)
        readset = reads_by_kmer[kmer]
        for read1, read2 in itertools.combinations(readset, 2):
            read_graph.add_edge(read1, read2)
    return read_graph


def write_partitions(read_graph, reads, ccprefix, logstream):
    """Given a read graph, write distinct partitions to separate files."""
    n = 0
    reads_in_ccs = 0
    cclog = open(ccprefix + '.cc.log', 'w')
    ccs = sorted(networkx.connected_components(read_graph), reverse=True,
                 # Sort first by number of reads, then by read names
                 key=lambda c: (len(c), sorted(c)))
    for n, cc in enumerate(ccs):
        print('CC', n, len(cc), cc, sep='\t', file=cclog)
        reads_in_ccs += len(cc)
        outfilename = '{:s}.cc{:d}.augfastq.gz'.format(ccprefix, n)
        with kevlar.open(outfilename, 'w') as outfile:
            for readid in cc:
                record = reads[readid]
                kevlar.print_augmented_fastx(record, outfile)
    message = '[kevlar::overlap] grouped {:d} reads'.format(reads_in_ccs)
    message += ' into {:d} connected components'.format(n + 1)
    print(message, file=logstream)
