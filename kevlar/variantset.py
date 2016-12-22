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
from sys import stdout
import kevlar


def setpair():
    return (set(), set())


class VariantSet(object):
    def __init__(self):
        self.linear_paths = defaultdict(setpair)
        self.kmers = defaultdict(setpair)
        self.read_ids = defaultdict(setpair)

    def add_kmer(self, kmer, read_id, linear_path):
        minkmer = kevlar.revcommin(kmer)
        minlp = kevlar.revcommin(linear_path)

        self.kmers[minkmer][0].add(minlp)
        self.kmers[minkmer][1].add(read_id)

        self.read_ids[read_id][0].add(minlp)
        self.read_ids[read_id][1].add(minkmer)

        self.linear_paths[minlp][0].add(read_id)
        self.linear_paths[minlp][1].add(minkmer)

    def collapse(self):
        upaths = dict()
        for path in sorted(self.linear_paths, key=len, reverse=True):
            pathrc = kevlar.revcom(path)
            merge = False
            for upath in upaths:
                if path in upath or pathrc in upath:
                    mergedreads = self.linear_paths[upath][0].union(
                        self.linear_paths[path][0]
                    )
                    mergedkmers = self.linear_paths[upath][1].union(
                        self.linear_paths[path][1]
                    )
                    self.linear_paths[upath] = (mergedreads, mergedkmers)
                    del self.linear_paths[path]
                    merge = True
                    break
            if merge is False:
                upaths[path] = True

    @property
    def npaths(self):
        return len(self.linear_paths)

    @property
    def nkmers(self):
        return len(self.kmers)

    @property
    def nreads(self):
        return len(self.read_ids)

    def iter_kmer(self):
        for minkmer in sorted(self.kmers):
            maxkmer = kevlar.revcom(minkmer)
            pathset = self.kmers[minkmer][0]
            readset = self.kmers[minkmer][1]
            yield minkmer, maxkmer, pathset, readset

    def iter_path(self):
        for minpath in sorted(self.linear_paths):
            maxpath = kevlar.revcom(minpath)
            readset = self.linear_paths[minpath][0]
            kmerset = self.linear_paths[minpath][1]
            yield minpath, maxpath, readset, kmerset

    def iter_read(self):
        for read_id in sorted(self.read_ids):
            pathset = self.read_ids[read_id][0]
            kmerset = self.read_ids[read_id][1]
            yield read_id, pathset, kmerset

    def kmer_table(self, outstream=stdout):
        for minkmer, maxkmer, pathset, readset in self.iter_kmer():
            kmerstr = ','.join((minkmer, maxkmer))
            pathstr = ','.join(sorted(pathset))
            readstr = ','.join(sorted(readset))
            print(kmerstr, pathstr, readstr, sep='\t', file=outstream)

    def path_table(self, outstream=stdout):
        for minpath, maxpath, readset, kmerset in self.iter_path():
            pathstr = ','.join((minpath, maxpath))
            readstr = ','.join(sorted(readset))
            kmerstr = ','.join(sorted(kmerset))
            print(pathstr, readstr, kmerstr, sep='\t', file=outstream)

    def read_table(self, outstream=stdout):
        for readid, pathset, kmerset in self.iter_read():
            pathstr = ','.join(sorted(pathset))
            kmerstr = ','.join(sorted(kmerset))
            print(readid, pathstr, kmerstr, sep='\t', file=outstream)
