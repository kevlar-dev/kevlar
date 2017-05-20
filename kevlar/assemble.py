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
        if logstream and n % 10000 == 0:  # pragma: no cover
            print('[kevlar::assemble]    loaded {:d} reads'.format(n),
                  file=logstream)
        reads[record.name] = record
        for kmer in record.ikmers:
            kmerseq = kevlar.revcommin(kmer.sequence)
            kmers[kmerseq].add(record.name)
    return reads, kmers


def graph_init_abund_check_pass(numreads, minkmer, minabund=5, maxabund=500,
                                logstream=None):
    """Check whether the k-mer falls within the expected range of abundance."""
    if maxabund and numreads > maxabund:
        msg = '            skipping k-mer with abundance {:d}'.format(numreads)
        print(msg, file=logstream)
        return False
    if minabund and numreads < minabund:
        message = '[kevlar::assemble] WARNING: k-mer {}'.format(minkmer)
        message += ' (rev. comp. {})'.format(kevlar.revcom(minkmer))
        message += ' only has abundance {:d}'.format(len(readnames))
        out = logstream if logstream is not None else sys.stderr
        print(message, file=out)
        return False
    return True


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
                       orient=pair.sameorient, tail=tailname)


def graph_init(reads, kmers, minabund=5, maxabund=500, logstream=None):
    """
    Initialize the "shared interesting k-mer" read graph.

    Iterate through each interesting k-mer, consider every pair of reads
    containing that k-mer, and determine whether there should be a edge between
    the pair of reads in the graph.
    """
    graph = networkx.Graph()
    nkmers = len(kmers)
    for n, minkmer in enumerate(kmers, 1):
        if n % 100 == 0:  # pragma: no cover
            msg = 'processed {:d}/{:d} shared novel k-mers'.format(n, nkmers)
            print('[kevlar::assemble]    ', msg, sep='', file=logstream)

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


def fetch_largest_overlapping_pair(graph, reads):
    """
    Grab the edge with the largest overlap in the graph.

    Sort the edges using 3 criteria. The first is the primary criterion, the
    other two ensure deterministic behavior.
        - overlap (largest first)
        - lexicographically smaller read name
        - lexicographically larger read name
    """
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
    return kevlar.overlap.OverlappingReadPair(
        tail=reads[read1],
        head=reads[read2],
        offset=graph[read1][read2]['offset'],
        overlap=graph[read1][read2]['overlap'],
        sameorient=graph[read1][read2]['orient'],
    )


def assemble_with_greed(reads, kmers, graph, debugout=None):
    """Find shortest common superstring using a greedy assembly algorithm."""
    count = 0
    while len(graph.edges()) > 0:
        count += 1

        pair = fetch_largest_overlapping_pair(graph, reads)
        newname = 'contig{:d}'.format(count)
        newrecord = merge_and_reannotate(pair, newname)
        if debugout:
            print('### DEBUG', pair.tail.name, pair.head.name, pair.offset,
                  pair.overlap, pair.sameorient, file=debugout)
            kevlar.print_augmented_fastq(newrecord, debugout)
        for kmer in newrecord.ikmers:
            kmerseq = kevlar.revcommin(kmer.sequence)
            for readname in kmers[kmerseq]:
                already_merged = readname not in graph
                current_contig = readname in [
                    pair.tail.name, pair.head.name, newname
                ]
                if already_merged or current_contig:
                    continue
                otherrecord = reads[readname]
                newpair = kevlar.overlap.calc_offset(
                    newrecord, otherrecord, kmerseq, debugout
                )
                if newpair == kevlar.overlap.INCOMPATIBLE_PAIR:
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
        graph.remove_node(pair.tail.name)
        graph.remove_node(pair.head.name)


def main(args):
    debugout = None
    if args.debug:
        debugout = args.logfile

    reads, kmers = load_reads(kevlar.open(args.augfastq, 'r'), debugout)
    inputreads = list(reads)
    graph = graph_init(reads, kmers, 2, args.max_abund, debugout)
    if args.gml:
        networkx.write_gml(graph, args.gml)
        message = '[kevlar::assemble] graph written to {}'.format(args.gml)
        print(message, file=args.logfile)

    ccs = list(networkx.connected_components(graph))
    assert len(ccs) == 1

    assemble_with_greed(reads, kmers, graph, debugout)

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
