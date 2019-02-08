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
from gzip import open as gzopen
from os import makedirs
from os.path import dirname
import re
from subprocess import Popen, PIPE, check_call
import sys
from tempfile import TemporaryFile

# Third-party libraries
import khmer
import networkx
import pysam

# Internal modules
from kevlar import seqio
from kevlar import sketch
from kevlar import reference
from kevlar import cigar
from kevlar import varmap
from kevlar import vcf
from kevlar import evaluate
from kevlar.intervalforest import IntervalForest
from kevlar.readpair import ReadPair
from kevlar.readgraph import ReadGraph
from kevlar.mutablestring import MutableString
from kevlar.sequence import parse_augmented_fastx, print_augmented_fastx
from kevlar.sequence import revcom
from kevlar.seqio import parse_partitioned_reads, parse_single_partition
from kevlar.timer import Timer
from kevlar.progress import ProgressIndicator

# Subcommands and command-line interface
from kevlar import novel
from kevlar import filter
from kevlar import augment
from kevlar import mutate
from kevlar import assemble
from kevlar import count
from kevlar import partition
from kevlar import localize
from kevlar import call
from kevlar import alac
from kevlar import varfilter
from kevlar import simlike
from kevlar import dist
from kevlar import split
from kevlar import unband
from kevlar import gentrio
from kevlar import cli

# C extensions
from kevlar.alignment import contig_align as align
import kevlar.assembly
import kevlar.sequence

from kevlar._version import get_versions
__version__ = get_versions()['version']
del get_versions


logstream = None
teelog = False


def plog(*args, **kwargs):
    """Print logging output."""
    if logstream is not None:
        print(*args, **kwargs, file=logstream)
    if logstream is None or teelog:
        print(*args, **kwargs, file=sys.stderr)


def open(filename, mode):
    if mode not in ('r', 'w'):
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
    networkx.write_gml(graph, outfilename, stringizer=str)
    message = '[kevlar] graph written to {}'.format(outfilename)
    print(message, file=logfile)


def multi_file_iter_khmer(filenames):
    for filename in filenames:
        for record in khmer.ReadParser(filename):
            yield record


def parse_bed(instream):
    for line in instream:
        if line.startswith('#'):
            continue
        line = line.strip()
        if line == '':
            continue
        values = re.split(r'\s+', line)
        chrom, start, end, *data = values
        yield chrom, int(start), int(end), data


def bedstream(bedfilelist):
    for bedfile in bedfilelist:
        fh = kevlar.open(bedfile, 'r')
        for values in parse_bed(fh):
            yield values


def vcf_header(outstream, version='4.2', source='kevlar', infoheader=False):
    print('##fileformat=VCFv', version, sep='', file=outstream)
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
