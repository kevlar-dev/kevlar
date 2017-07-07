#!/usr/bin/env python
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
import pandas
import networkx
import screed
import khmer
import kevlar
from kevlar.seqio import load_reads_and_kmers


def merge_pair(pair):
    """
    Assemble a pair of overlapping reads.

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


def fetch_largest_overlapping_pair(graph, reads):
    """
    Grab the edge with the largest overlap in the graph.

    Sort the edges using 3 criteria. The first is the primary criterion, the
    other two ensure deterministic behavior.
        - FIXME
        - overlap (largest first)
        - lexicographically smaller read name
        - lexicographically larger read name
    """
    edges = sorted(
        graph.edges(),
        reverse=True,
        key=lambda e: (
            sum(graph.degree([e[0], e[1]]).values()),
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
        swapped=graph[read1][read2]['swapped'],
    )


def assemble_with_greed(reads, kmers, graph, ccindex, debugout=None):
    """Find shortest common superstring using a greedy assembly algorithm."""
    count = 0
    while len(graph.edges()) > 0:
        count += 1

        pair = fetch_largest_overlapping_pair(graph, reads)
        newname = 'contig{:d};cc={:d}'.format(count, ccindex)
        newrecord = merge_and_reannotate(pair, newname)
        if debugout:
            print('### DEBUG', pair.tail.name, pair.head.name, pair.offset,
                  pair.overlap, pair.sameorient, file=debugout)
            kevlar.print_augmented_fastx(newrecord, debugout)
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
                                   orient=newpair.sameorient, tail=tn,
                                   swapped=newpair.swapped)
            kmers[kmerseq].add(newrecord.name)
        reads[newrecord.name] = newrecord
        graph.add_node(newrecord.name)
        graph.remove_node(pair.tail.name)
        graph.remove_node(pair.head.name)


def prune_graph(graph, quant=0.1):
    edge_adj_deg = list()
    for edge in graph.edges_iter():
        degree_dict = graph.degree([edge[0], edge[1]])
        agg_degree = sum(degree_dict.values())
        edge_adj_deg.append(agg_degree)

    edges = pandas.DataFrame(edge_adj_deg)
    threshold = edges[0].quantile(quant)
    edges_to_drop = list()  # Don't remove edges while iterating through them
    for edge in graph.edges_iter():
        degree_dict = graph.degree([edge[0], edge[1]])
        agg_degree = sum(degree_dict.values())
        if agg_degree < threshold:
            edges_to_drop.append(edge)

    for edge in edges_to_drop:
        graph.remove_edge(edge[0], edge[1])

    return len(edges_to_drop)


def main(args):
    debugout = None
    if args.debug:
        debugout = args.logfile

    reads, kmers = load_reads_and_kmers(kevlar.open(args.augfastq, 'r'),
                                        args.logfile)
    inputreads = list(reads)
    message = 'loaded {:d} reads'.format(len(inputreads))
    message += ' and {:d} interesting k-mers'.format(len(kmers))
    print('[kevlar::assemble]', message, file=args.logfile)

    graph = kevlar.overlap.graph_init_strict(reads, kmers, args.min_abund,
                                             args.max_abund, debugout)
    message = 'initialized "shared interesting k-mers" graph'
    message += ' with {:d} nodes'.format(graph.number_of_nodes())
    message += ' and {:d} edges'.format(graph.number_of_edges())
    # If number of nodes is less than number of reads, it's probably because
    # some reads have no valid overlaps with other reads.
    print('[kevlar::assemble]', message, file=args.logfile)

    if args.gml:
        tempgraph = graph.copy()
        for n1, n2 in tempgraph.edges():
            ikmerset = tempgraph[n1][n2]['ikmers']
            ikmerstr = ','.join(ikmerset)
            tempgraph[n1][n2]['ikmers'] = ikmerstr
        networkx.write_gml(tempgraph, args.gml)
        message = '[kevlar::assemble] graph written to {}'.format(args.gml)
        print(message, file=args.logfile)

    edges_dropped = prune_graph(graph)
    cc_stream = networkx.connected_component_subgraphs(graph)
    ccs = [cc for cc in cc_stream if cc.number_of_edges() > 0]
    ccnodes = sum([cc.number_of_nodes() for cc in ccs])
    message = 'dropped {:d} edges'.format(edges_dropped)
    message += ', graph now has {:d} connected component(s)'.format(len(ccs))
    message += ', {:d} nodes'.format(ccnodes)
    message += ', and {:d} edges'.format(graph.number_of_edges())
    print('[kevlar::assemble]', message, file=args.logfile)
    if len(ccs) > 1:
        message = 'multiple connected components designated by cc=N in output'
        print('[kevlar::assemble] WARNING:', message, file=args.logfile)

    contigcount = 0
    unassembledcount = 0
    outstream = kevlar.open(args.out, 'w')
    for n, cc in enumerate(ccs, 1):
        assemble_with_greed(reads, kmers, cc, n, debugout)
        for seqname in cc.nodes():
            if seqname in inputreads:
                unassembledcount += 1
                continue
            contigcount += 1
            contigrecord = reads[seqname]
            kevlar.print_augmented_fastx(contigrecord, outstream)

    assembledcount = ccnodes - unassembledcount
    message = '[kevlar::assemble] assembled'
    message += ' {:d}/{:d} reads'.format(assembledcount, ccnodes)
    message += ' from {:d} connected component(s)'.format(len(ccs))
    message += ' into {:d} contig(s)'.format(contigcount)
    print(message, file=args.logfile)
