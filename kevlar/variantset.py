#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
import kevlar
import screed


def setpair():
    return (set(), set())


def revcommin(seq):
    revcom = screed.dna.reverse_complement(str(seq))
    minseq = sorted((seq, revcom))[0]
    return minseq


class VariantSet(obj):
    def __init__(self):
        self.linear_paths = defaultdict(setpair)
        self.kmers = defaultdict(setpair)
        self.read_ids = defaultdict(setpair)

    def add_kmer(kmer, read_id, linear_path):
        minkmer = revcommin(kmer)
        minlp = revcommin(linear_path)

        self.kmers[minkmer][0].add(minlp)
        self.kmers[minkmer][1].add(read_id)

        self.read_ids[read_id][0].add(minlp)
        self.read_ids[read_id][1].add(minkmer)

        self.linear_paths[minlp][0].add(readid)
        self.linear_paths[minlp][1].add(minkmer)

    def iter_kmer():
        for minkmer in self.kmers:
            maxkmer = screed.dna.reverse_complement(minkmer)
            pathset = self.kmers[minkmer][0]
            readset = self.kmers[minkmer][1]
            yield minkmer, maxkmer, pathset, readset

    def iter_path():
        for minpath in self.linear_paths:
            maxpath = screed.dna.reverse_complement(minpath)
            readset = self.kmers[minpath][0]
            kmerset = self.kmers[minpath][1]
            yield minpath, maxpath, readset, kmerset

    def iter_read():
        for read_id in self.read_ids:
            pathset = self.read_ids[read_id][0]
            kmerset = self.read_ids[read_id][1]
            return readid, pathset, kmerset
