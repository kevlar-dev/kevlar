#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

try:
    import __builtin__ as builtins
except:  # pragma: no cover
    import builtins
from collections import namedtuple
from kevlar import seqio
from kevlar.seqio import parse_augmented_fastq, print_augmented_fastq
from kevlar import dump
from kevlar import find
from kevlar import collect
from kevlar import filter
from kevlar import reaugment
from kevlar import assemble
from kevlar import cli
from kevlar.variantset import VariantSet
from kevlar.timer import Timer
from gzip import open as gzopen
import khmer
import screed

from kevlar._version import get_versions
__version__ = get_versions()['version']
del get_versions


def open(filename, mode):
    if mode not in ['r', 'w']:
        raise ValueError('invalid mode "{}"'.format(mode))
    openfunc = builtins.open
    if filename.endswith('.gz'):
        openfunc = gzopen
        mode += 't'
    return openfunc(filename, mode)


def calc_fpr(table):
    """Stolen shamelessly from khmer/__init__.py"""
    sizes = table.hashsizes()
    n_ht = float(len(sizes))
    occupancy = float(table.n_occupied())
    min_size = min(sizes)
    fp_one = occupancy / min_size
    fp_all = fp_one ** n_ht
    return fp_all


def revcom(seq):
    return screed.dna.reverse_complement(str(seq))


def revcommin(seq):
    rc = revcom(seq)
    minseq = sorted((seq, rc))[0]
    return minseq


def same_seq(seq1, seq2, seq2revcom=None):
    if seq2revcom is None:
        seq2revcom = revcom(seq2)
    return seq1 == seq2 or seq1 == seq2revcom


def load_sketch(filename, count=False, graph=False, smallcount=False):
    if count and graph:
        if smallcount:
            createfunc = khmer._SmallCountgraph
        else:
            createfunc = khmer._Countgraph
    elif count and not graph:
        if smallcount:
            createfunc = khmer._SmallCounttable
        else:
            createfunc = khmer._Counttable
    elif not count and graph:
        createfunc = khmer._Nodegraph
    elif not count and not graph:
        createfunc = khmer._Nodetable

    sketch = createfunc(1, [1])
    sketch.load(filename)
    return sketch


KmerOfInterest = namedtuple('KmerOfInterest', 'sequence offset abund')
