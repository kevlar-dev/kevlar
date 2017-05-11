#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict, namedtuple
import argparse
import itertools
import sys
try:
    import networkx
except:
    pass
import screed
import khmer
import kevlar


OverlappingReadPair = namedtuple('OverlappingReadPair',
                                 'tail head offset overlap sameorient')
IncompatiblePair = OverlappingReadPair(None, None, None, None, None)


def subparser(subparsers):
    subparser = subparsers.add_parser('assemble')
    subparser.add_argument('-d', '--debug', action='store_true',
                           help='print debugging output')
    subparser.add_argument('-o', '--out', metavar='FILE',
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('--gml', metavar='FILE',
                           help='write graph to .gml file')
    subparser.add_argument('-x', '--max-abund', type=int, metavar='X',
                           default=500, help='discard interesting k-mers that '
                           'occur more than X times')
    subparser.add_argument('augfastq', help='annotated reads in augmented '
                           'Fastq format')


def print_read_pair(read1, pos1, read2, pos2, ksize, offset, overlap,
                    sameorient, outstream):
    seq2 = read2.sequence if sameorient else kevlar.revcom(read2.sequence)
    print('\nDEBUG shared interesting kmer: ', read1.name,
          ' --(overlap={:d})'.format(overlap),
          '(offset={:d})({})--> '.format(offset, sameorient), read2.name, '\n',
          read1.sequence, '\n', ' ' * pos1, '|' * ksize, '\n', ' ' * offset,
          seq2, '\n', sep='', file=outstream)


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
    maxkmer = kevlar.revcom(minkmer)
    kmer1 = [k for k in read1.ikmers
             if kevlar.same_seq(k.sequence, minkmer, maxkmer)][0]
    kmer2 = [k for k in read2.ikmers
             if kevlar.same_seq(k.sequence, minkmer, maxkmer)][0]
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

    headseq = head.sequence if sameorient else kevlar.revcom(head.sequence)
    seg2offset = len(head.sequence) - len(tail.sequence) + offset
    if offset + len(headseq) <= len(tail.sequence):
        segment1 = tail.sequence[offset:offset+len(headseq)]
        segment2 = headseq
        seg2offset = None
    else:
        segment1 = tail.sequence[offset:]
        segment2 = headseq[:-seg2offset]

    overlap1 = len(segment1)
    overlap2 = len(segment2)
    if overlap1 != overlap2:
        print(
            'DEBUG '
            'tail="{tail}" head="{head}" offset={offset} altoffset={altoffset}'
            ' tailoverlap={overlap} headoverlap={headover} tailolvp={tailseq}'
            ' headolvp={headseq} kmer={minkmer},{maxkmer} tailseq={tailread}'
            ' headseq={headread}'.format(
                tail=read1.name, head=read2.name, offset=offset,
                altoffset=seg2offset, overlap=overlap1,
                headover=len(segment2), tailseq=segment1, headseq=segment2,
                minkmer=minkmer, maxkmer=maxkmer, tailread=tail.sequence,
                headread=head.sequence
            ), file=sys.stderr
        )
    assert overlap1 == overlap2
    if segment1 != segment2:
        return IncompatiblePair

    if debugstream:
        print_read_pair(tail, tailpos, head, headpos, ksize, offset, overlap1,
                        sameorient, debugstream)

    return OverlappingReadPair(tail=tail, head=head, offset=offset,
                               overlap=overlap1, sameorient=sameorient)


def merge_pair(pair):
    """
    Assemble a pair of overlapping reads.

    Given a pair of compatible overlapping reads, collapse and merge them into
    a single sequence.
    """
    headseq = pair.head.sequence
    if pair.sameorient is False:
        headseq = kevlar.revcom(pair.head.sequence)
    headindex = len(pair.tail.sequence) - pair.offset
    headsuffix = headseq[headindex:]
    tailprefix = pair.tail.sequence[pair.offset:pair.offset+pair.overlap]
    assert tailprefix == headseq[:headindex], \
        'error: attempted to assemble incompatible reads'

    if headseq in pair.tail.sequence:
        return pair.tail.sequence
    return pair.tail.sequence + headsuffix


def merge_and_reannotate(pair, newname):
    """
    Assemble a pair of overlapping reads and resolve their interesting k-mers.

    When a pair of compatible reads is merged, the offset of the interesting
    k-mers must be computed for one of the reads.
    """
    contig = merge_pair(pair)
    newrecord = screed.Record(name=newname, sequence=contig,
                              ikmers=pair.tail.ikmers)
    ksize = len(pair.tail.ikmers[0].sequence)
    if pair.sameorient:
        minoffset2keep = len(pair.tail.sequence) - pair.offset - ksize
        keepers = [ik for ik in pair.head.ikmers if ik.offset > minoffset2keep]
        for k in keepers:
            ikmer = kevlar.KmerOfInterest(k.sequence, k.offset + pair.offset,
                                          k.abund)
            newrecord.ikmers.append(ikmer)
    else:
        maxoffset2keep = pair.offset - ksize
        keepers = [ik for ik in pair.head.ikmers if ik.offset < maxoffset2keep]
        for k in keepers:
            ikmer = kevlar.KmerOfInterest(
                kevlar.revcom(k.sequence),
                len(pair.head.sequence) - k.offset - ksize + pair.offset,
                k.abund,
            )
            newrecord.ikmers.append(ikmer)

    return newrecord


def load_reads(instream, logstream=None):
    """
    Load reads into lookup tables for convenient access.

    The first table is a dictionary of reads indexed by read name, and the
    second table is a dictionary of read sets indexed by an interesting k-mer.
    """
    reads = dict()
    kmers = defaultdict(set)
    for n, record in enumerate(kevlar.parse_augmented_fastq(instream), 1):
        if logstream and n % 10000 == 0:
            print('[kevlar::assemble]    loaded {:d} reads'.format(n),
                  file=logstream)
        reads[record.name] = record
        for kmer in record.ikmers:
            kmerseq = kevlar.revcommin(kmer.sequence)
            kmers[kmerseq].add(record.name)
    return reads, kmers


def graph_init(reads, kmers, maxabund=500, logstream=None):
    graph = networkx.Graph()
    nkmers = len(kmers)
    for n, minkmer in enumerate(kmers, 1):
        if n % 100 == 0:
            msg = 'processed {:d}/{:d} shared novel k-mers'.format(n, nkmers)
            print('[kevlar::assemble]    ', msg, sep='', file=logstream)
        readnames = kmers[minkmer]
        if maxabund and len(readnames) > maxabund:
            msg = '            skipping k-mer with abundance {:d}'.format(
                len(readnames)
            )
            print(msg, file=logstream)
            continue
        if len(readnames) < 2:
            message = '[kevlar::assemble] WARNING: k-mer {}'.format(minkmer)
            message += ' (rev. comp. {})'.format(kevlar.revcom(minkmer))
            message += ' only has abundance {:d}'.format(len(readnames))
            out = logstream if logstream is not None else sys.stderr
            print(message, file=out)
            continue
        readset = [reads[rn] for rn in readnames]
        for read1, read2 in itertools.combinations(readset, 2):
            pair = calc_offset(read1, read2, minkmer, logstream)
            if pair is IncompatiblePair:  # Shared k-mer but bad overlap
                continue
            tailname, headname = pair.tail.name, pair.head.name
            if tailname in graph and headname in graph[tailname]:
                assert graph[tailname][headname]['offset'] == pair.offset
                if graph[tailname][headname]['tail'] == tailname:
                    assert graph[tailname][headname]['overlap'] == pair.overlap
            else:
                graph.add_edge(tailname, headname, offset=pair.offset,
                               overlap=pair.overlap, ikmer=minkmer,
                               orient=pair.sameorient, tail=tailname)
    return graph


def assemble_with_greed(reads, kmers, graph):
    pass


def main(args):
    debugout = None
    if args.debug:
        debugout = args.logfile

    reads, kmers = load_reads(kevlar.open(args.augfastq, 'r'), debugout)
    inputreads = list(reads.keys())
    graph = graph_init(reads, kmers, args.max_abund, debugout)
    if args.gml:
        networkx.write_gml(graph, args.gml)
        message = '[kevlar::assemble] graph written to {}'.format(args.gml)
        print(message, file=args.logfile)

    ccs = list(networkx.connected_components(graph))
    assert len(ccs) == 1

    count = 0
    while len(graph.edges()) > 0:
        count += 1
        # Sort the edges using 3 criteria. The first is the primary criterion,
        # the other two ensure deterministic behavior.
        #     - overlap (largest first)
        #     - lexicographically smaller read name
        #     - lexicographically larger read name
        edges = sorted(
            graph.edges(),
            reverse=True,
            key=lambda e: (
                graph[e[0]][e[1]]['overlap'],
                max(e),
                min(e),
            )
        )
        read1, read2 = edges[0]  # biggest overlap (greedy algorithm)
        if read2 == graph[read1][read2]['tail']:
            read1, read2 = read2, read1
        pair = OverlappingReadPair(
            tail=reads[read1],
            head=reads[read2],
            offset=graph[read1][read2]['offset'],
            overlap=graph[read1][read2]['overlap'],
            sameorient=graph[read1][read2]['orient'],
        )
        newname = 'contig{:d}'.format(count)
        newrecord = merge_and_reannotate(pair, newname)
        if debugout:
            print('### DEBUG', read1, read2, pair.offset, pair.overlap,
                  pair.sameorient, file=debugout)
            kevlar.print_augmented_fastq(newrecord, debugout)
        for kmer in newrecord.ikmers:
            kmerseq = kevlar.revcommin(kmer.sequence)
            for readname in kmers[kmerseq]:
                already_merged = readname not in graph
                current_contig = readname in [read1, read2, newname]
                if already_merged or current_contig:
                    continue
                otherrecord = reads[readname]
                newpair = calc_offset(
                    newrecord, otherrecord, kmerseq, debugout
                )
                if newpair == IncompatiblePair:
                    continue
                tn, hn = newpair.tail.name, newpair.head.name
                if tn in graph and hn in graph[tn]:
                    assert graph[tn][hn]['overlap'] == newpair.overlap
                    if graph[tn][hn]['tail'] == newpair.tail:
                        assert graph[tn][hn]['offset'] == newpair.offset
                else:
                    graph.add_edge(tn, hn, offset=newpair.offset,
                                   overlap=newpair.overlap, ikmer=kmerseq,
                                   orient=newpair.sameorient, tail=tn)
            kmers[kmerseq].add(newrecord.name)
        reads[newrecord.name] = newrecord
        graph.add_node(newrecord.name)
        graph.remove_node(read1)
        graph.remove_node(read2)

    contigcount = 0
    unassembledcount = 0
    outstream = kevlar.open(args.out, 'w')
    for seqname in graph.nodes():
        if seqname in inputreads:
            unassembledcount += 1
            continue
        contigcount += 1
        contigrecord = reads[seqname]
        kevlar.print_augmented_fastq(contigrecord, outstream)

    assembledcount = len(inputreads) - unassembledcount
    message = '[kevlar::assemble] assembled'
    message += ' {:d}/{:d} reads'.format(assembledcount, len(inputreads))
    message += ' into {:d} contig(s)'.format(contigcount)
    print(message, file=args.logfile)
