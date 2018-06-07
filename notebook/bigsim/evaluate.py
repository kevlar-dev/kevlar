#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
from collections import defaultdict
import sys

import intervaltree
from intervaltree import IntervalTree
import kevlar
from kevlar.vcf import VCFReader


class IntervalForest(object):
    """Single point of access for a labeled set of interval trees.

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
    index = IntervalForest()
    with kevlar.open(filename, 'r') as instream:
        for line in instream:
            if line.strip() == '':
                continue
            values = line.strip().split()
            position = int(values[0])
            if values[1] == 'Del':
                position -= int(values[2])
            strrepr = '{}<-{}'.format(values[1], values[2])
            index.insert(chrlabel, position, position + 1, strrepr)
    return index


def assess_variants(variants, index, delta=10):
    variants_by_class = defaultdict(list)
    correct = set()
    false = set()
    missing = set()
    found = set()
    mapping = defaultdict(list)

    for varcall in variants:
        if varcall.filterstr != 'PASS':
            continue
        callclass = varcall.attribute('CALLCLASS')
        if callclass is not None:
            variants_by_class[callclass].append(varcall)
            continue
        hits = index.query(varcall.seqid, varcall.position, delta=delta)
        if hits == set():
            false.add(varcall)
        else:
            correct.add(varcall)
            found.update(hits)
            for hit in hits:
                mapping[hit].append(varcall)

    for callclass, calllist in variants_by_class.items():
        nmatches = 0
        match = None
        localfound = set()
        for varcall in calllist:
            hits = index.query(varcall.seqid, varcall.position, delta=delta)
            if hits == set():
                continue
            else:
                nmatches += 1
                match = varcall
                localfound = hits
        if nmatches == 0:
            false.add(calllist[0])
        else:
            assert nmatches > 0, nmatches
            if nmatches > 1:
                print('WARNING: found', nmatches, 'matches for CALLCLASS',
                      callclass, file=sys.stderr)
            correct.add(varcall)
            found.update(localfound)
            for hit in localfound:
                mapping[hit].append(varcall)

    for label, interval in index:
        if interval not in found:
            missing.add(interval)

    return correct, false, missing, mapping


def parse_mvf(instream, skip=1):
    if skip:
        [next(instream) for _ in range(skip)]
    for line in instream:
        values = line.split('\t')
        seqid = values[0]
        position = int(values[1])
        yield seqid, position, line


def assess_variants_mvf(variants, index, delta=10):
    variants_by_class = defaultdict(list)
    correct = set()
    false = set()
    missing = set()
    found = set()
    mapping = defaultdict(list)

    for seqid, position, line in variants:
        hits = index.query(seqid, position, delta=delta)
        if hits == set():
            false.add(line)
        else:
            correct.add(line)
            found.update(hits)
            for hit in hits:
                mapping[hit].append(line)

    for label, interval in index:
        if interval not in found:
            missing.add(interval)

    return correct, false, missing, mapping


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tolerance', type=int, metavar='T', default=10,
                        help='extend real variants by T nucleotides when '
                        'querying for overlap with variant calls; default is '
                        '10')
    parser.add_argument('--mvf', action='store_true', help='input is in .mvf '
                        'format; default is .vcf')
    parser.add_argument('--correct', help='print correct variants to file')
    parser.add_argument('--missing', help='print missing variants to file')
    parser.add_argument('--false', help='print false variants to file')
    parser.add_argument('--collisions', help='print calls that match the same '
                        'variant')
    parser.add_argument('simvar', help='simulated variants (in custome 3-4 '
                        'column tabular format)')
    parser.add_argument('varcalls', help='VCF file of variant calls')
    args = parser.parse_args()

    index = populate_index_from_simulation(args.simvar, 'chr17')
    if args.mvf:
        reader = parse_mvf(kevlar.open(args.varcalls, 'r'))
        assess_func = assess_variants_mvf
    else:
        reader = VCFReader(kevlar.open(args.varcalls, 'r'))
        assess_func = assess_variants
    correct, false, missing, mapping = assess_func(
        reader, index, delta=args.tolerance
    )

    numcollisions = 0
    for variant, calllist in mapping.items():
        if len(calllist) > 1:
            numcollisions += 1
    if numcollisions > 0:
        print('WARNING:', numcollisions, 'variants matched by multiple calls',
              file=sys.stderr)
        if args.collisions:
            with open(args.collisions, 'w') as outstream:
                for variant, calllist in mapping.items():
                    if len(calllist) > 1:
                        print('\n#VARIANT:', variant, file=outstream)
                        for varcall in calllist:
                            if args.mvf:
                                print('    -', varcall, end='', file=outstream)
                            else:
                                print('    -', varcall.vcf, file=outstream)

    if args.missing:
        with open(args.missing, 'w') as outstream:
            for variant in missing:
                print(variant.begin, *variant.data.split('<-'), sep='\t',
                      file=outstream)

    if args.correct:
        outstream = kevlar.open(args.correct, 'w')
        if args.mvf:
            for varcall in correct:
                print(varcall, end='', file=outstream)
        else:
            writer = kevlar.vcf.VCFWriter(outstream)
            for varcall in correct:
                writer.write(varcall)

    if args.false:
        outstream = kevlar.open(args.false, 'w')
        if args.mvf:
            for varcall in false:
                print(varcall, end='', file=outstream)
        else:
            writer = kevlar.vcf.VCFWriter(outstream)
            for varcall in false:
                writer.write(varcall)

    print('Correct:', len(mapping))
    print('False:', len(false))
    print('Missing:', len(missing))
