#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict
import itertools
import networkx
import kevlar


class ReadGraph(networkx.Graph):

    def __init__(self, data=None, **attr):
        """
        Constructor

        The only addition to the base class is a dictionary to store sets of
        reads containing each "interesting" k-mer. These k-mers are used to
        build out the graph edges.
        """
        self.ikmers = defaultdict(set)
        super(ReadGraph, self).__init__(data, **attr)

    def load(self, instream, minabund=None, maxabund=None, dedup=False):
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

        for record in kevlar.parse_augmented_fastx(instream):
            minread = kevlar.revcommin(record.sequence)
            if minread in unique_reads:
                continue
            unique_reads.add(minread)

            self.add_node(record.name, record=record)
            for kmer in record.ikmers:
                kmerseq = kevlar.revcommin(kmer.sequence)
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
                          orient=pair.sameorient, tail=tailname,
                          swapped=pair.swapped)

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
                    record1 = self[read1]['record']
                    record2 = self[read2]['record']
                    pair = kevlar.overlap.calc_offset(record1, record2, kmer)
                    if pair is kevlar.overlap.INCOMPATIBLE_PAIR:
                        # Shared k-mer but bad overlap
                        continue
                    self.check_edge(pair, kmer)
                else:
                    self.add_edge(read1, read2)
