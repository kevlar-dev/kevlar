#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict
from Bio.trie import trie
from sys import stdout
import kevlar


class VariantSet(object):
    def __init__(self):
        self.contigs = defaultdict(set)  # contig -> set(kmers)
        self.kmers = defaultdict(set)  # kmer -> set(read IDs)
        self._kmer_instances = 0
        self._reads = 0

    def add_kmer(self, kmer, read_id):
        min_kmer = kevlar.revcommin(kmer)
        self.kmers[min_kmer].add(read_id)
        self._kmer_instances += 1

    def add_contig(self, contig, kmer):
        min_contig = kevlar.revcommin(contig)
        min_kmer = kevlar.revcommin(kmer)
        self.contigs[min_contig].add(min_kmer)

    def collapse(self):
        unique_contigs = list()
        suffix_tree = trie()
        for contig in sorted(self.contigs, key=len, reverse=True):
            contigrc = kevlar.revcom(contig)
            query = suffix_tree.with_prefix(contig) + \
                suffix_tree.with_prefix(contigrc)
            if query == []:
                unique_contigs.append(contig)
                for i in range(len(contig)):
                    suffix_tree[contig[i:]] = 1
            else:
                longest = sorted(query, key=len, reverse=True)[0]
                mergedkmers = self.contigs[longest].union(
                    self.contigs[contig]
                )
                self.contigs[longest] = mergedkmers
                del self.contigs[contig]

    @property
    def ncontigs(self):
        return len(self.contigs)

    @property
    def nkmers(self):
        return len(self.kmers)

    @property
    def nkmerinst(self):
        return self._kmer_instances

    @property
    def nreads(self):
        reads = set()
        for kmer in self.kmers:
            kreads = self.kmers[kmer]
            reads = reads.union(kreads)
        return len(reads)

    def __iter__(self):
        for mincontig in sorted(self.contigs):
            maxcontig = kevlar.revcom(mincontig)
            kmers = self.contigs[mincontig]
            reads = set()
            for kmer in kmers:
                reads = reads.union(self.kmers[kmer])
            yield mincontig, maxcontig, kmers, reads

    def write(self, outstream=stdout):
        print('Contig,ContigRevCom\tNumReads\tNumKmers\tReads\tKmers',
              file=outstream)
        for mincontig, maxcontig, kmerset, readset in self:
            contigstr = ','.join((mincontig, maxcontig))
            readstr = ','.join(sorted(readset))
            kmerstr = ','.join(sorted(kmerset))
            print(contigstr, len(readset), len(kmerset), readstr, kmerstr,
                  sep='\t', file=outstream)
