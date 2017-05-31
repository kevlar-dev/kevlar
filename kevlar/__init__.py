#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
try:
    import __builtin__ as builtins
except:  # pragma: no cover
    import builtins
from collections import namedtuple
import re
import sys
from kevlar import seqio
from kevlar import overlap
from kevlar.seqio import parse_augmented_fastq, print_augmented_fastq
from kevlar import dump
from kevlar import novel
from kevlar import collect
from kevlar import filter
from kevlar import reaugment
from kevlar import mutate
from kevlar import assemble
from kevlar import count
from kevlar import partition
from kevlar import counting
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
    if filename in ['-', None]:
        filehandle = sys.stdin if mode == 'r' else sys.stdout
        return filehandle
    openfunc = builtins.open
    if filename.endswith('.gz'):
        openfunc = gzopen
        mode += 't'
    return openfunc(filename, mode)


def sketch_autoload(infile, count=True, graph=False,
                    ksize=31, table_size=1e4, num_tables=4,
                    num_bands=0, band=0):
    """
    Use file extension to conditionally load sketch into memory.

    If the file extension is one of the following, treat the file as a sketch
    that has been written to disk and load it with `kevlar.load_sketch`.
    Sketch attributes such as ksize, table size, and number of tables will
    be set automatically.

    - `.ct` or `.counttable`: `Counttable`
    - `.nt` or `.nodetable`: `Nodetable`
    - `.cg` or `.countgraph`: `Countgraph`
    - `.ng` or `.nodegraph`: `Nodegraph`

    Otherwise, a sketch will be created using the specified arguments and the
    input file will be treated as a Fasta/Fastq file to be loaded with
    `.consume_seqfile` or `.consume_seqfile_banding`.
    """
    sketch_extensions = (
        '.ct', '.counttable', '.nt', '.nodetable',
        '.cg', '.countgraph', '.ng', '.nodegraph',
    )

    if infile.endswith(sketch_extensions):
        return load_sketch(infile, count=count, graph=graph, smallcount=False)
    else:
        sketch = allocate_sketch(ksize, table_size, num_tables, count=count,
                                 graph=graph, smallcount=False)
        if num_bands > 1:
            sketch.consume_seqfile_banding(infile, num_bands, band)
        else:
            sketch.consume_seqfile(infile)
        return sketch


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
    """Convenience function for loading a sketch from the specified file."""
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


def allocate_sketch(ksize, target_tablesize, num_tables=4, count=False,
                    graph=False, smallcount=False):
    """Convenience function for allocating memory for a new sketch."""
    if count and graph:
        if smallcount:
            createfunc = khmer.SmallCountgraph
        else:
            createfunc = khmer.Countgraph
    elif count and not graph:
        if smallcount:
            createfunc = khmer.SmallCounttable
        else:
            createfunc = khmer.Counttable
    elif not count and graph:
        createfunc = khmer.Nodegraph
    elif not count and not graph:
        createfunc = khmer.Nodetable

    sketch = createfunc(ksize, target_tablesize, num_tables)
    return sketch


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
