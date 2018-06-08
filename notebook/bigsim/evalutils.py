#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

"""Please destroy this software after publication. kthxbye.

This code was written to facilitate evaluations for a scientific paper and is
not intended to be general purpose. It is shared in the spirit of transparency,
not as exemplary code. For more thoughts on this, see the post and comments at
http://ivory.idyll.org/blog/2015-how-should-we-think-about-research-software.html.
"""

import argparse
from collections import defaultdict
import sys

import intervaltree
from intervaltree import IntervalTree
import kevlar
from kevlar.vcf import VCFReader


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


def populate_index_from_simulation(filename, chrlabel):
    """Index loading function

    The file containing simulated variants is a tab-separated file with 3 or 4
    values in each row.

        498149  Ins     208
        498767  Del     343
        677694  Del     9
        691785  t       a       SNV

    The first column is the position of the variant (on chr17 in this case).
    For indels, the variant type is in the second column and the length is in
    the third column. For SNVs, the alternate allele is in the second column,
    the reference allele is in the third column, and the variant type is in the
    fourth column.

    For some reason, deletion variants are listed by their terminal nucleotide
    instead of by their initial nucleotide, thus the correction marked by the
    comment "compute first nucleotide".
    """
    index = IntervalForest()
    with kevlar.open(filename, 'r') as instream:
        for line in instream:
            if line.strip() == '':
                continue
            values = line.strip().split()
            position = int(values[0])
            if values[1] == 'Del':  # compute first nucleotide
                position -= int(values[2])
            strrepr = '{}<-{}'.format(values[1], values[2])
            index.insert(chrlabel, position, position + 1, strrepr)
    return index


def compact(reader, index, delta=10):
    """Compact variants by call class

    Variant calls labeled with the same `CALLCLASS` attribute were predicted
    from the same partition (set of reads). While more than one of these may
    indeed be true variants with respect to the reference, we expect only one
    *de novo* variant per partition.

    This function assumes the variant calls are sorted by likelihood score (the
    `LIKESCORE` attribute in kevlar). Any call that does not pass filters is
    ignored. Then, for each CALLCLASS with multiple passing calls, all calls
    are discarded except for the one matching the true variant. If the
    CALLCLASS has no calls matching a true variant, all of the calls are
    discarded except for the highest scoring call.
    """
    variants_by_class = defaultdict(list)
    calls = list()
    for varcall in reader:
        if varcall.filterstr != 'PASS':
            continue
        callclass = varcall.attribute('CALLCLASS')
        if callclass is None:
            calls.append(varcall)
        else:
            variants_by_class[callclass].append(varcall)

    for callclass, calllist in variants_by_class.items():
        nmatches = 0
        match = None
        for varcall in calllist:
            hits = index.query(varcall.seqid, varcall.position, delta=delta)
            if hits == set():
                continue
            else:
                nmatches += 1
                if match is None:
                    match = varcall
                    localfound = hits
        if nmatches == 0:
            calls.append(calllist[0])
        else:
            assert nmatches > 0, nmatches
            if nmatches > 1:
                print('WARNING: found', nmatches, 'matches for CALLCLASS',
                      callclass, file=sys.stderr)
            calls.append(match)

    calls.sort(key=lambda c: float(c.attribute('LIKESCORE')), reverse=True)
    return calls



def assess_variants(variants, index, delta=10):
    """Assess kevlar calls in .vcf format."""
    variants_by_class = defaultdict(list)
    correct = set()
    false = set()
    missing = set()
    found = set()
    mapping = defaultdict(list)

    for varcall in variants:
        assert varcall.filterstr == 'PASS'
        hits = index.query(varcall.seqid, varcall.position, delta=delta)
        if hits == set():
            false.add(varcall)
        else:
            correct.add(varcall)
            found.update(hits)
            assert len(hits) == 1
            for hit in hits:
                mapping[hit].append(varcall)

    for label, interval in index:
        if interval not in found:
            missing.add(interval)

    return correct, false, missing, mapping


def assess_variants_mvf(variants, index, delta=10):
    """Assess GATK calls in .mvf format."""
    variants_by_class = defaultdict(list)
    correct = set()
    false = set()
    missing = set()
    found = set()
    mapping = defaultdict(list)

    for rowindex, row in variants.iterrows():
        hits = index.query(row['CHROM'], row['POS'], delta=delta)
        varstr = '{:s}:{:d}({:s})'.format(row['CHROM'], row['POS'], row['CHILD_GT'])
        if hits == set():
            false.add(varstr)
        else:
            correct.add(varstr)
            found.update(hits)
            assert len(hits) == 1
            for hit in hits:
                mapping[hit].append(varstr)

    for label, interval in index:
        if interval not in found:
            missing.add(interval)

    return correct, false, missing, mapping
