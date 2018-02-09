#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

# Core libraries
import builtins
from collections import namedtuple
from gzip import open as gzopen
from os import makedirs
from os.path import dirname
import re
import sys

# Third-party libraries
import khmer
import screed

# Internal modules
from kevlar import seqio
from kevlar import overlap
from kevlar import sketch
from kevlar import vcf
from kevlar.mutablestring import MutableString
from kevlar.readgraph import ReadGraph
from kevlar.seqio import parse_augmented_fastx, print_augmented_fastx
from kevlar.seqio import parse_partitioned_reads, parse_single_partition
from kevlar.variantset import VariantSet
from kevlar.timer import Timer

# Subcommands and command-line interface
from kevlar import dump
from kevlar import novel
from kevlar import filter
from kevlar import augment
from kevlar import mutate
from kevlar import assemble
from kevlar import count
from kevlar import effcount
from kevlar import partition
from kevlar import localize
from kevlar import call
from kevlar import alac
from kevlar import simplex
from kevlar import simlike
from kevlar import split
from kevlar import gentrio
from kevlar import cli

# C extensions
from kevlar.alignment import contig_align as align
import kevlar.assembly

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


def mkdirp(path, trim=False):
    outdir = dirname(path) if trim else path
    makedirs(outdir, exist_ok=True)
    return outdir


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


def vcf_header(outstream, version='4.2', source='kevlar', infoheader=False):
    print('##fileFormat=VCFv', version, sep='', file=outstream)
    print('##source=', source, sep='', file=outstream)
    if infoheader:
        print('##INFO=<GT,Number=3,Type=String,Description="Genotypes of each '
              'individual in the trio (proband, mother, father)">',
              file=outstream)
    print('##INFO=<VW,Number=1,Type=String,Description="Genomic interval '
          'bounding all k-mers that contain the alternate allele">',
          file=outstream)
    print('##INFO=<RW,Number=1,Type=String,Description="Genomic interval '
          'bounding all k-mers that contain the reference allele">',
          file=outstream)
    print('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
          sep='\t', file=outstream)


KmerOfInterest = namedtuple('KmerOfInterest', 'sequence offset abund')
