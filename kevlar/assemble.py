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


SequenceOverlap = namedtuple('SequenceOverlap', 'seq1 seq2 offset sameorient')


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
    if pos2 > pos1:
        tail, head = read2, read1
        offset = pos2 - pos1

    segment1 = tail.sequence[offset:]
    headseq = head.sequence if sameorient else kevlar.revcom(head.sequence)
    segment2 = headseq[:(len(headseq)-offset)]
    if segment1 != segment2:
        return None, None, None, None

    if debugstream:
        print('\nDEBUG shared interesting kmer: ', tail.name,
              ' --({:d})({})--> '.format(offset, sameorient), head.name, '\n',
              tail.sequence, '\n', ' ' * pos1, '|' * ksize, '\n', ' ' * offset,
              headseq, '\n', sep='', file=debugstream)

    return tail, head, offset, sameorient


def collapse(seq1, seq2, offset, sameorient):
    if sameorient is False:
        seq2 = kevlar.revcom(seq2)
    firstnewnucl = len(seq2) - (len(seq2) - len(seq1) + offset)
    seq2suffix = seq2[firstnewnucl:]
    return seq1 + seq2suffix


def annotate_contig(record1, record2, newcontig, newname, offset, sameorient):
    assert offset > 1
    newrecord = screed.Record(name=newname, sequence=newcontig,
                              ikmers=record1.ikmers)
    ksize = len(record1.ikmers[0].sequence)
    if sameorient:
        # Next two lines are at high risk for off-by-one errors
        minoffset2keep = len(record2.sequence) - offset - ksize
        keepers = [ik for ik in record2.ikmers if ik.offset > minoffset2keep]
        for k in keepers:
            ikmer = kevlar.KmerOfInterest(k.sequence, k.offset + offset,
                                          k.abund)
            newrecord.ikmers.append(ikmer)
    else:
        maxoffset2keep = offset + ksize
        keepers = [ik for ik in record2.ikmers if ik.offset < maxoffset2keep]
        for k in keepers:
            print()
            ikmer = kevlar.KmerOfInterest(
                kevlar.revcom(k.sequence),
                len(record2.sequence) - k.offset - ksize + offset,
                k.abund,
            )
            newrecord.ikmers.append(ikmer)

    return newrecord


def main(args):
    reads = dict()            # key: read ID, value: record
    kmers = defaultdict(set)  # key: k-mer (min repr), value: set of read IDs
    for n, record in enumerate(kevlar.parse_augmented_fastq(args.augfastq), 1):
        if n % 10000 == 0:
            print('[kevlar::assemble]    loaded {:d} reads'.format(n),
                  file=args.logfile)
        reads[record.name] = record
        for kmer in record.ikmers:
            kmerseq = kevlar.revcommin(kmer.sequence)
            kmers[kmerseq].add(record.name)

    debugout = None
    if args.debug:
        debugout = args.logfile
    graph = networkx.Graph()  # DiGraph?
    nkmers = len(kmers)
    for n, minkmer in enumerate(kmers, 1):
        if n % 100 == 0:
            msg = 'processed {:d}/{:d} shared novel k-mers'.format(n, nkmers)
            print('[kevlar::assemble]    ', msg, sep='', file=args.logfile)
        readnames = kmers[minkmer]
        if len(readnames) > args.max_abund:
            msg = '            skipping k-mer with abundance {:d}'.format(
                len(readnames)
            )
            print(msg, file=sys.stderr)
            continue
        assert len(readnames) > 1
        readset = [reads[rn] for rn in readnames]
        for read1, read2 in itertools.combinations(readset, 2):
            tail, head, offset, sameorient = calc_offset(read1, read2, minkmer,
                                                         debugout)
            if tail is None:  # Shared k-mer but bad overlap
                continue
            if tail.name in graph and head.name in graph[tail.name]:
                assert graph[tail.name][head.name]['offset'] == offset
            else:
                graph.add_edge(tail.name, head.name, offset=offset,
                               ikmer=minkmer, orient=sameorient)

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

            newcontig = collapse(seq1, seq2, offset, sameorient)
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
