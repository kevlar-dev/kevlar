#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

# Core libraries
from __future__ import print_function
try:
    import __builtin__ as builtins
except:
    import builtins
from collections import namedtuple
from gzip import open as gzopen
import re
import sys

# Third-party libraries
import khmer
import screed

# Internal modules
from kevlar import seqio
from kevlar import overlap
from kevlar import counting
from kevlar import sketch
from kevlar.seqio import parse_augmented_fastx, print_augmented_fastx
from kevlar.variantset import VariantSet
from kevlar.timer import Timer

# Subcommands and command-line interface
from kevlar import dump
from kevlar import novel
from kevlar import collect
from kevlar import filter
from kevlar import reaugment
from kevlar import mutate
from kevlar import assemble
from kevlar import count
from kevlar import partition
from kevlar import localize
from kevlar import call
from kevlar import cli

# C extension(s)
from kevlar.alignment import contig_align as align

from kevlar._version import get_versions
__version__ = get_versions()['version']
del get_versions


def open(filename, mode):
    if mode not in ['r', 'w']:
        raise ValueError('invalid mode "{}"'.format(mode))
    if filename in ['-', None]:
        filehandle = sys.stdin if mode == 'r' else sys.stdout
        return filehandle
    openfunc = builtins.open
    if filename.endswith('.gz'):
        openfunc = gzopen
        mode += 't'
    return openfunc(filename, mode)


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


def to_gml(graph, outfilename, logfile=sys.stderr):
    """Write the given read graph to a GML file."""
    if not outfilename.endswith('.gml'):
        print('[kevlar] WARNING: GML files usually need extension .gml',
              file=logfile)
    networkx.write_gml(graph, outfilename)
    message = '[kevlar] graph written to {}'.format(args.gml)
    print(message, file=logfile)


def multi_file_iter_screed(filenames):
    for filename in filenames:
        for record in screed.open(filename):
            yield record


def multi_file_iter_khmer(filenames):
    for filename in filenames:
        for record in khmer.ReadParser(filename):
            yield record


def clean_subseqs(sequence, ksize):
    for subseq in re.split('[^ACGT]', sequence):
        if len(subseq) >= ksize:
            yield subseq


KmerOfInterest = namedtuple('KmerOfInterest', 'sequence offset abund')
