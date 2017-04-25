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
import networkx
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
    subparser.add_argument('-o', '--out', metavar='FILE', default=sys.stdout,
                           type=argparse.FileType('w'),
                           help='output file; default is terminal (stdout)')
    subparser.add_argument('--gml', metavar='FILE',
                           help='write graph to .gml file')
    subparser.add_argument('-x', '--max-abund', type=int, metavar='X',
                           default=500, help='discard interesting k-mers that '
                           'occur more than X times')
    subparser.add_argument('augfastq', type=argparse.FileType('r'),
                           help='annotated reads in augmented Fastq format')


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
    offset = pos1 - pos2
    if pos2 > pos1 or len(read2.sequence) > len(read1.sequence):
        tail, head = read2, read1
        offset = pos2 - pos1

    segment1 = tail.sequence[offset:]
    bpoverlap = len(segment1)
    headseq = head.sequence if sameorient else kevlar.revcom(head.sequence)
    segment2 = headseq[:(len(headseq)-offset)]
    assert len(segment2) == bpoverlap
    if segment1 != segment2:
        return IncompatiblePair

    if debugstream:
        print('\nDEBUG shared interesting kmer: ', tail.name,
              ' --({:d})({})--> '.format(offset, sameorient), head.name, '\n',
              tail.sequence, '\n', ' ' * pos1, '|' * ksize, '\n', ' ' * offset,
              headseq, '\n', sep='', file=debugstream)

    return OverlappingReadPair(tail=tail, head=head, offset=offset,
                               overlap=bpoverlap, sameorient=sameorient)


def merge_pair(pair):
    """
    Assemble a pair of overlapping reads.

    Given a pair of compatible overlapping reads, collapse and merge them into
    a single sequence.
    """
    assert len(pair.tail.sequence) >= len(pair.head.sequence)
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
            print()
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
            print('[kevlar::assemble]    ', msg, sep='', file=args.logfile)
        readnames = kmers[minkmer]
        if maxabund and len(readnames) > maxabund:
            msg = '            skipping k-mer with abundance {:d}'.format(
                len(readnames)
            )
            print(msg, file=logstream)
            continue
        assert len(readnames) > 1
        readset = [reads[rn] for rn in readnames]
        for read1, read2 in itertools.combinations(readset, 2):
            pair = calc_offset(read1, read2, minkmer, logstream)
            if pair is IncompatiblePair:  # Shared k-mer but bad overlap
                continue
            tailname, headname = pair.tail.name, pair.head.name
            if tailname in graph and headname in graph[tailname]:
                assert graph[tailname][headname]['offset'] == pair.offset
            else:
                graph.add_edge(tailname, headname, offset=pair.offset,
                               overlap=pair.overlap, ikmer=minkmer,
                               orient=pair.sameorient)
    return graph


def main(args):
    debugout = None
    if args.debug:
        debugout = args.logfile

    reads, kmers = load_reads(args.augfastq, debugout)
    graph = graph_init(reads, kmers, args.max_abund, debugout)

    reads_assembled = 0
    ccs = list(networkx.connected_components(graph))
    for n, cc in enumerate(ccs, 1):
        count = 0
        while len(graph.edges()) > 0:
            count += 1
            edges = sorted(graph.edges(),
                           key=lambda e: graph[e[0]][e[1]]['offset'])
            read1, read2 = edges[0]  # biggest overlap
            seq1 = reads[read1].sequence
            seq2 = reads[read2].sequence
            offset = graph[read1][read2]['offset']
            sameorient = graph[read1][read2]['orient']

            if offset == 0:
                graph.remove_node(read2)
                continue

            newcontig = merge(seq1, seq2, offset, sameorient)
            newname = 'contig{:d}'.format(count)
            newrecord = annotate_contig(reads[read1], reads[read2], newcontig,
                                        newname, offset, sameorient)
            if debugout:
                print('# DEBUG', read1, read2, offset, sameorient,
                      file=debugout)
                kevlar.print_augmented_fastq(newrecord, debugout)

            for kmer in newrecord.ikmers:
                kmerseq = kevlar.revcommin(kmer.sequence)
                for readname in kmers[kmerseq]:
                    if readname not in graph:
                        continue
                    otherrecord = reads[readname]
                    tail, head, poffset, sameorient = calc_offset(
                        newrecord, otherrecord, kmerseq, debugout
                    )
                    if tail is None:
                        continue
                    if tail.name in graph and head.name in graph[tail.name]:
                        assert graph[tail.name][head.name]['offset'] == poffset
                    else:
                        graph.add_edge(tail.name, head.name, offset=poffset,
                                       ikmer=minkmer, orient=sameorient)

            kmers[kmerseq].add(newrecord.name)
            graph.remove_node(read1)
            graph.remove_node(read2)

        reads_assembled += len(cc)
        print('CC', n, len(cc), sep='\t', file=args.out)

    if args.gml:
        networkx.write_gml(graph, args.gml)
        message = '[kevlar::assemble] graph written to {}'.format(args.gml)
        print(message, file=args.logfile)

    message = '[kevlar::assemble] {:d} reads'.format(reads_assembled)
    message += ' assembled into {:d} contigs'.format(n)
    print(message, file=args.logfile)
