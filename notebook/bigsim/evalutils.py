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
ivory.idyll.org/blog/2015-how-should-we-think-about-research-software.html.
"""

import argparse
from collections import defaultdict
import sys

import intervaltree
from intervaltree import IntervalTree
import kevlar
from kevlar.vcf import VCFReader
import pandas


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
        if end is None:
            end = start
        if delta > 0:
            start -= delta
            end += delta
        if end == start:
            return self.trees[label][start]
        else:
            return self.trees[label][start:end]


def populate_index_from_simulation(instream, chrlabel):
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


def populate_index_from_bed(instream):
    index = IntervalForest()
    for line in instream:
        if line.strip() == '':
            continue
        values = line.strip().split()
        chrom = values[0]
        start, end = [int(coord) for coord in values[1:3]]
        if start == end:
            start -= 1
        index.insert(chrom, start, end)
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
    calls = [c for c in calls if float(c.attribute('LIKESCORE')) > 0.0]
    return calls


def load_kevlar_vcf(filename, index, delta=10, vartype=None, minlength=None, maxlength=None):
    reader = kevlar.vcf.VCFReader(kevlar.open(filename, 'r'))
    reader.suppress_filter_warnings = True
    calls = compact(reader, index, delta=delta)
    if vartype:
        calls = subset_vcf(calls, vartype, minlength=minlength, maxlength=maxlength)
    return calls


def load_scalpel_vcf(filename, vartype=None, minlength=None, maxlength=None, cov='30'):
    reader = kevlar.vcf.VCFReader(kevlar.open(filename, 'r'))
    reader.suppress_filter_warnings = True
    calls = list(reader)
    calls.sort(key=lambda c: float(c.attribute('CHI2')), reverse=True)
    if vartype:
        calls = subset_vcf(calls, vartype, minlength=minlength, maxlength=maxlength)
    return calls


def load_triodenovo_vcf(filename, vartype=None, minlength=None, maxlength=None, cov='30'):
    reader = kevlar.vcf.VCFReader(kevlar.open(filename, 'r'))
    reader.suppress_filter_warnings = True
    calls = list(reader)
    calls.sort(key=lambda c: float(c._sample_data['Kid_'+ cov +'x']['DQ']), reverse=True)
    if vartype:
        calls = subset_vcf(calls, vartype, minlength=minlength, maxlength=maxlength)
    return calls


def load_gatk_mvf(filename, vartype=None, minlength=None, maxlength=None):
    calls = pandas.read_table(filename, sep='\t').sort_values('CHILD_DP', ascending=False)
    if vartype:
        calls = subset_mvf(calls, vartype, minlength=minlength, maxlength=maxlength)
    return calls


def assess_variants_vcf(variants, index, delta=10):
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
        varstr = '{:s}:{:d}({:s})'.format(
            row['CHROM'], row['POS'], row['CHILD_GT']
        )
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


def subset_variants(variants, vartype, minlength=None, maxlength=None):
    assert vartype in ('SNV', 'INDEL')
    for line in variants:
        if line.strip() == '':
            continue
        values = line.strip().split()
        vtype = 'SNV' if values[-1] == 'SNV' else 'INDEL'
        if vtype != vartype:
            continue

        if vartype == 'SNV':
            yield line
            continue

        indellength = int(values[2])
        if minlength and indellength < minlength:
            continue
        if maxlength and indellength > maxlength:
            continue
        yield line


def subset_variants_bed(variants, vartype, minlength=None, maxlength=None):
    assert vartype in ('SNV', 'INDEL')
    for line in variants:
        if line.strip() == '':
            continue
        values = line.strip().split()
        if len(values) == 5:
            if vartype == 'SNV':
                yield line
            continue
        if vartype == 'INDEL':
            if len(values) > 3:
                indellength = int(values[3])
            else:
                indellength = int(values[2]) - int(values[1])
            if minlength and indellength < minlength:
                continue
            if maxlength and indellength > maxlength:
                continue
            yield line


def subset_vcf(varcalls, vartype, minlength=None, maxlength=None):
    assert vartype in ('SNV', 'INDEL')
    for call in varcalls:
        refrlen = len(call._refr)
        altlen = len(call._alt)
        calltype = 'SNV'
        if refrlen > 1 or altlen > 1:
            calltype = 'INDEL'
        if calltype != vartype:
            continue

        if vartype == 'SNV':
            yield call
            continue

        indellength = abs(refrlen - altlen)
        if minlength and indellength < minlength:
            continue
        if maxlength and indellength > maxlength:
            continue
        yield call


def subset_mvf(varcalls, vartype, minlength=None, maxlength=None):
    assert vartype in ('SNV', 'INDEL')
    def determine_type(genotype):
        g1, g2 = genotype.replace('|', '/').split('/')
        if len(g1) == 1 and len(g2) == 1:
            return 'SNV'
        return 'INDEL'

    def determine_length(genotype):
        g1, g2 = genotype.replace('|', '/').split('/')
        if len(g1) == 1 and len(g2) == 1:
            return 0
        return abs(len(g1) - len(g2))

    varcalls['VARTYPE'] = varcalls['CHILD_GT'].apply(determine_type)
    varcalls['VARLENGTH'] = varcalls['CHILD_GT'].apply(determine_length)

    varcalls = varcalls[varcalls['VARTYPE'] == vartype]
    if minlength:
        varcalls = varcalls[varcalls['VARLENGTH'] >= minlength]
    if maxlength:
        varcalls = varcalls[varcalls['VARLENGTH'] <= maxlength]
    return varcalls
