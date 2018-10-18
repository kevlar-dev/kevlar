#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
import itertools
import networkx
import kevlar


class ReadGraph(networkx.Graph):

    def __init__(self, data=None, **attr):
        """
        Constructor

        In addition to the base class, we add a dictionary to store sets of
        reads containing each "interesting" k-mer. These k-mers are used to
        build out the graph edges.

        Also, we store the names of the input reads so that reads with no
        connections to other reads can be distinguished from assembled contigs.
        """
        self.ikmers = defaultdict(set)
        self.readnames = set()
        super(ReadGraph, self).__init__(data, **attr)

    def full_cc(self, cc):
        sg = self.subgraph(cc).copy()
        sg = ReadGraph(data=sg)
        sg.ikmers = self.ikmers
        sg.readnames = self.readnames
        return sg

    def get_record(self, recordname):
        return self.node[recordname]['record']

    def load(self, readstream, minabund=None, maxabund=None, dedup=False):
        """
        Load reads and interesting k-mers into a graph structure.

        A graph node is created for each read, and a set of reads containing
        each interesting k-mer is stored. If abundance thresholds are enforced,
        do a second in-memory pass over the k-mers to discard any that don't
        satisfy the threshold criteria.

        Set `dedup=True` to deduplicate read sequences for handling PCR
        duplicates. This doesn't do a proper check (i.e. check both pairs
        against genome), but simply makes sure that only one copy of each read
        sequence is loaded. This is implemented with a very naive and resource
        intensive approach, so this mode should only be used on small (e.g.
        already partitioned) graphs.
        """
        temp_ikmers = defaultdict(set)
        unique_reads = set()

        for record in readstream:
            if dedup:
                minread = kevlar.revcommin(record.sequence)
                if minread in unique_reads:
                    continue
                unique_reads.add(minread)

            self.add_node(record.name, record=record)
            self.readnames.add(record.name)
            for kmer in record.annotations:
                kmerseq = kevlar.revcommin(record.ikmerseq(kmer))
                temp_ikmers[kmerseq].add(record.name)

        if minabund is None and maxabund is None:
            self.ikmers = temp_ikmers
        else:
            for kmer in temp_ikmers:
                readset = temp_ikmers[kmer]
                abund = len(readset)
                minfail = minabund and abund < minabund
                maxfail = maxabund and abund > maxabund
                if not minfail and not maxfail:
                    self.ikmers[kmer] = readset

    def check_edge(self, pair, minkmer):
        """
        Add edge between 2 nodes in the "shared interesting k-mer" read graph.

        If the edge already exists, make sure that the existing edge matches
        the edge that would have been added.
        """
        tailname, headname = pair.tail.name, pair.head.name
        if tailname in self and headname in self[tailname]:
            assert self[tailname][headname]['offset'] == pair.offset
            if self[tailname][headname]['tail'] == tailname:
                assert self[tailname][headname]['overlap'] == pair.overlap
            self[tailname][headname]['ikmers'].add(minkmer)
        else:
            self.add_edge(tailname, headname, offset=pair.offset,
                          overlap=pair.overlap, ikmers=set([minkmer]),
                          orient=pair.sameorient, tail=tailname)

    def populate_edges(self, strict=False):
        """
        Instantiate edges between nodes that share an interesting k-mer.

        Setting `strict=True` will result in edges between reads only when they
        have a perfect match in their overlap.
        """
        for kmer in self.ikmers:
            readset = self.ikmers[kmer]
            for read1, read2 in itertools.combinations(readset, 2):
                if strict:
                    if read1 in self and read2 in self[read1]:
                        continue
                    record1 = self.get_record(read1)
                    record2 = self.get_record(read2)
                    pair = kevlar.ReadPair(record1, record2, kmer)
                    if pair.incompatible:
                        # Shared k-mer but bad overlap
                        continue
                    self.check_edge(pair, kmer)
                else:
                    self.add_edge(read1, read2)

    def partitions(self, dedup=True, minabund=None, maxabund=None,
                   abundfilt=False):
        """
        Retrieve all partitions (connected components) from this graph.

        The `minabund` and `maxabund` parameters are used at graph construction
        time to filter out k-mers whose abundance is too large or small. If
        `abundfilt` is true, the minimum bundance is also applied to the number
        of sequences (reads or contigs) in the partition.
        """
        for cc in sorted(networkx.connected_components(self), reverse=True,
                         # Sort first by number of reads, then by read names
                         key=lambda c: (len(c), sorted(c))):
            if len(cc) == 1 and list(cc)[0] in self.readnames:
                continue  # Skip unassembled input reads
            if dedup:
                partition = ReadGraph()
                readstream = [self.get_record(readid) for readid in cc]
                partition.load(readstream, minabund, maxabund, dedup=True)
                assert partition.number_of_nodes() > 0
                if abundfilt:
                    if minabund and partition.number_of_nodes() < minabund:
                        continue  # Skip partitions that are too small
                # # Ill-advised strategy for pre-emptively discarding
                # # unassemblable partitions. It turns out that number of edges
                # # in the read graph (even in strict mode) doesn't really
                # # distinguish between what assembles and what doesn't.
                # if partition.number_of_nodes() < 10:
                #     partition.populate_edges(strict=True)
                #     nedges = partition.number_of_edges()
                #     if minabund and nedges < minabund:
                #         continue
                yield partition
            else:
                yield cc
