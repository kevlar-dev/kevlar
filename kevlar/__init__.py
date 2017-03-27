#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from . import fasta
from . import dump
from . import find
from . import collect
from .variantset import VariantSet
from .timer import Timer
import screed

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


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


def parse_augmented_fastq(instream):
    record = None
    annot_kmers = dict()

    for line in instream:
        if line.startswith('@'):
            if record is not None:
                yield record, annot_kmers
                annot_kmers = dict()
            readid = line[1:].strip()
            seq = next(instream).strip()
            _ = next(instream)
            qual = next(instream).strip()
            record = screed.Record(name=readid, sequence=seq, quality=qual)
        elif line.endswith('#\n'):
            offset = len(line) - len(line.lstrip())
            line = line.strip()[:-1]
            abundances = re.split('\s+', line)
            kmer = abundances.pop(0)
            annot_kmers[offset] = (kmer, abundances)
    if record is not None:
        yield record, annot_kmers
