#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019 Battelle National Biodefense Institute
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
from intervaltree import IntervalTree


class IntervalForest(object):
    """Single point of access for a labeled set of interval trees.

    By default, queries return only exact matches (delta=0). For inexact
    matches, set delta=X for X nucleotide extensions in each direction
    for the query.

    >>> index = IntervalForest()
    >>> index.insert('chr17', 238026, 238046)
    >>> index.insert('chr17', 1533596, 1533597, 'C->A')
    >>> index.query('chr17', 1533500, 1533700)
    {Interval(1533596, 1533597, 'C->A')}
    >>> index.query('chr17', 238006)
    set()
    >>> index.query('chr17', 238006, delta=30)
    {Interval(238026, 238046, 'chr17:238026-238046')}
    >>> index.query('chr4', 1533500, 1533700)
    set()
    """

    def __init__(self):
        self.trees = defaultdict(IntervalTree)

    def __iter__(self):
        for label, tree in sorted(self.trees.items()):
            for interval in tree:
                yield label, interval

    def __len__(self):
        return sum([len(tree) for tree in self.trees.values()])

    def __iter__(self):
        for label, tree in self.trees.items():
            for interval in tree:
                yield interval.data

    def insert(self, label, start, end, data=None):
        assert label is not None
        if data is None:
            data = '{:s}:{:d}-{:d}'.format(label, start, end)
        self.trees[label][start:end] = data

    def query(self, label, start, end=None, delta=0):
        if label not in self.trees:
            return set()
        if delta > 0:
            if end:
                end += delta
            else:
                end = start + delta
            start -= delta
        if end is None:
            return self.trees[label][start]
        else:
            return self.trees[label][start:end]
