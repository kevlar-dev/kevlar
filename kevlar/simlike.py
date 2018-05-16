#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import khmer
import kevlar
from kevlar.vcf import Variant
from math import log
import scipy.stats


def get_abundances(sequence, case, controls, refr):
    """Create a nested list of k-mer abundances.

    abunds = [
        [15, 14, 13, 16, 14, 15, 14, 14],  # k-mer abundances from case/proband
        [0, 0, 1, 0, 2, 10, 0, 0],  # k-mer abundances from parent/control 1
        [0, 1, 1, 0, 1, 0, 2, 0],  # k-mer abundances from parent/control 2
    ]
    """
    kmers = case.get_kmers(sequence)
    valid_kmers = [k for k in kmers if refr.get(k) == 0]
    ndropped = len(kmers) - len(valid_kmers)

    abundances = list()
    for _ in range(len(controls) + 1):
        abundances.append(list())
    for kmer in valid_kmers:
        abund = case.get(kmer)
        abundances[0].append(abund)
        for control, abundlist in zip(controls, abundances[1:]):
            abund = control.get(kmer)
            abundlist.append(abund)
    return abundances, ndropped
