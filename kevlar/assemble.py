#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict, namedtuple
import itertools
import sys
import pandas
import networkx
import screed
import khmer
import kevlar
from kevlar.seqio import load_reads_and_kmers


# =============================================================================
# Fermi-lite assembly mode
# =============================================================================

def assemble_fml_asm(readstream, logstream=sys.stderr):
    reads = [r for r in readstream]
    assembler = kevlar.assembly.fml_asm(reads)
    for n, contig in enumerate(assembler, 1):
        name = 'contig{:d}'.format(n)
        record = screed.Record(name=name, sequence=contig)
        yield next(kevlar.augment.augment(reads, [record]))


def main_fml_asm(args):
    reads = kevlar.parse_augmented_fastx(kevlar.open(args.augfastq, 'r'))
    outstream = kevlar.open(args.out, 'w')
    for contig in assemble_fml_asm(reads):
        kevlar.print_augmented_fastx(contig, outstream)


# =============================================================================
# Greedy assembly mode
# =============================================================================


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


def fetch_largest_overlapping_pair(graph):
    """
    Grab the edge with the largest overlap in the graph.

    Sort the edges using 4 criteria. The first is the primary criterion, the
    other three ensure deterministic behavior.
        - the aggregate degree of the adjacent nodes
        - overlap (largest first)
        - lexicographically smaller read name
        - lexicographically larger read name
    """
    edges = sorted(
        graph.edges(),
        reverse=True,
        key=lambda e: (
            sum([d[1] for d in graph.degree([e[0], e[1]])]),
            graph[e[0]][e[1]]['overlap'],
            max(e),
            min(e),
        )
    )
    read1, read2 = edges[0]  # biggest overlap (greedy algorithm)
    if read2 == graph[read1][read2]['tail']:
        read1, read2 = read2, read1
    return kevlar.overlap.OverlappingReadPair(
        tail=graph.get_record(read1),
        head=graph.get_record(read2),
        offset=graph[read1][read2]['offset'],
        overlap=graph[read1][read2]['overlap'],
        sameorient=graph[read1][read2]['orient'],
        swapped=graph[read1][read2]['swapped'],
    )


def assemble_with_greed(graph, ccindex, debugout=None):
    """Find shortest common superstring using a greedy assembly algorithm."""
    count = 0
    while len(graph.edges()) > 0:
        count += 1

        pair = fetch_largest_overlapping_pair(graph)
        newname = 'contig{:d}:cc={:d}'.format(count, ccindex)
        newrecord = merge_and_reannotate(pair, newname)
        if debugout:  # pragma: no cover
            print('### DEBUG', pair.tail.name, pair.head.name, pair.offset,
                  pair.overlap, pair.sameorient, file=debugout)
            kevlar.print_augmented_fastx(newrecord, debugout)
        for kmer in newrecord.ikmers:
            kmerseq = kevlar.revcommin(kmer.sequence)
            for readname in graph.ikmers[kmerseq]:
                already_merged = readname not in graph
                current_contig = readname in [
                    pair.tail.name, pair.head.name, newname
                ]
                if already_merged or current_contig:
                    continue
                otherrecord = graph.get_record(readname)
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
            graph.ikmers[kmerseq].add(newrecord.name)
        graph.add_node(newrecord.name, record=newrecord)
        graph.remove_node(pair.tail.name)
        graph.remove_node(pair.head.name)


def prune_graph(graph, quant=0.1):
    edge_adj_deg = list()
    for edge in graph.edges:
        degree = graph.degree([edge[0], edge[1]])
        agg_degree = sum([d[1] for d in degree])
        edge_adj_deg.append(agg_degree)

    edges = pandas.DataFrame(edge_adj_deg)
    threshold = edges[0].quantile(quant)
    edges_to_drop = list()  # Don't remove edges while iterating through them
    for edge in graph.edges:
        degree = graph.degree([edge[0], edge[1]])
        agg_degree = sum(d[1] for d in degree)
        if agg_degree < threshold:
            edges_to_drop.append(edge)

    for edge in edges_to_drop:
        graph.remove_edge(edge[0], edge[1])

    return len(edges_to_drop)


def assemble_greedy(readstream, compat=0.6, debug=False, logstream=sys.stderr):
    debugout = None
    if debug:  # pragma: no cover
        debugout = logstream

    graph = kevlar.ReadGraph()
    graph.load(readstream)
    inputreads = set(graph.nodes())
    orignumnodes = graph.number_of_nodes()
    message = 'loaded {:d} reads'.format(orignumnodes)
    message += ' and {:d} interesting k-mers'.format(len(graph.ikmers))
    print('[kevlar::assemble::greedy]', message, file=logstream)

    graph.populate_edges(strict=True)
    message = 'populated "shared interesting k-mers" graph'
    message += ' with {:d} edges'.format(graph.number_of_edges())
    # If number of nodes is less than number of reads, it's probably because
    # some reads have no valid overlaps with other reads.
    print('[kevlar::assemble::greedy]', message, file=logstream)

    if graph.number_of_edges() == 0:
        print('[kevlar::assemble::greedy] nothing to be done', file=logstream)
        return

    edges_dropped = prune_graph(graph)
    cc_stream = networkx.connected_component_subgraphs(graph, copy=False)
    ccs = [cc for cc in cc_stream if cc.number_of_edges() > 0]
    ccnodes = sum([cc.number_of_nodes() for cc in ccs])
    message = 'dropped {:d} edges'.format(edges_dropped)
    message += ', graph now has {:d} connected component(s)'.format(len(ccs))
    message += ', {:d} nodes'.format(ccnodes)
    message += ', and {:d} edges'.format(graph.number_of_edges())
    print('[kevlar::assemble::greedy]', message, file=logstream)
    if ccnodes / orignumnodes < compat:
        msg = 'only {:d} of {:d} reads '.format(ccnodes, orignumnodes)
        msg += 'have compatible overlaps; discarding'
        print('[kevlar::assemble::greedy]', msg, file=logstream)
        return
    if len(ccs) > 1:
        message = 'multiple connected components designated by cc=N in output'
        print('[kevlar::assemble::greedy] WARNING:', message, file=logstream)

    contigs2report = list()
    unassembledcount = 0
    for n, cc in enumerate(ccs, 1):
        cc = graph.full_cc(cc)
        ccreads = list()
        for readname in cc.nodes():
            if readname in inputreads:
                ccreads.append(graph.get_record(readname))
        assemble_with_greed(cc, n, debugout)
        for seqname in cc.nodes():
            if seqname in inputreads:
                unassembledcount += 1
                continue
            contigrecord = cc.get_record(seqname)
            contig = next(kevlar.augment.augment(ccreads, [contigrecord]))
            contigs2report.append(contig)

    assembledcount = ccnodes - unassembledcount
    message = 'assembled {:d}/{:d} reads'.format(assembledcount, ccnodes)
    message += ' from {:d} connected component(s)'.format(len(ccs))
    message += ' into {:d} contig(s)'.format(len(contigs2report))
    print('[kevlar::assemble::greedy]', message, file=logstream)
    if assembledcount / ccnodes < compat:
        message = 'too few reads assembled; discarding'
        print('[kevlar::assemble::greedy]', message, file=logstream)
        return

    for contig in contigs2report:
        yield contig


def main_greedy(args):
    readstream = kevlar.parse_augmented_fastx(kevlar.open(args.augfastq, 'r'))
    outstream = kevlar.open(args.out, 'w')
    contigstream = assemble_greedy(readstream, debug=args.debug,
                                   logstream=args.logfile)
    for contig in contigstream:
        kevlar.print_augmented_fastx(contig, outstream)


# =============================================================================
# Main method
# =============================================================================

def main(args):
    mainfunc = main_fml_asm
    if args.greedy:
        mainfunc = main_greedy
    mainfunc(args)
