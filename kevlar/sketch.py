#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import khmer


sketch_constructor_by_filename_extension = {
    '.nt': khmer._Nodetable,
    '.ng': khmer._Nodegraph,
    '.ct': khmer._Counttable,
    '.cg': khmer._Countgraph,
    '.sct': khmer._SmallCounttable,
    '.scg': khmer._SmallCountgraph,
    '.nodetable': khmer._Nodetable,
    '.nodegraph': khmer._Nodegraph,
    '.counttable': khmer._Counttable,
    '.countgraph': khmer._Countgraph,
    '.smallcounttable': khmer._SmallCounttable,
    '.smallcountgraph': khmer._SmallCountgraph,
}


class KevlarSketchTypeError(ValueError):
    pass


def estimate_fpr(sketch):
    """
    Get a rough estimate of the false positive rate of this sketch.

    Stolen shamelessly from khmer/__init__.py
    """
    sizes = sketch.hashsizes()
    n_ht = float(len(sizes))
    occupancy = float(sketch.n_occupied())
    min_size = min(sizes)
    fp_one = occupancy / min_size
    fp_all = fp_one ** n_ht
    return fp_all


def load(filename):
    """
    Convenience function for loading a sketch from the specified file.

    Unfortunately, this relies on filename extensions, which are subject to
    human error. But until khmer stores all relevant information in the file
    itself and enables loading directly from file contents, this is the best we
    can do.
    """
    extensions = tuple(sketch_constructor_by_filename_extension)
    if not filename.endswith(extensions):
        message = 'unable to determine sketch type from filename ' + filename
        raise KevlarSketchTypeError(message)

    ext = '.' + filename.split('.')[-1]
    constructor = sketch_constructor_by_filename_extension[ext]

    sketch = constructor(1, [1])
    sketch.load(filename)
    return sketch


def allocate(ksize, target_tablesize, num_tables=4, count=False, graph=False,
             smallcount=False):
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
    else:
        assert not count and not graph
        createfunc = khmer.Nodetable

    sketch = createfunc(ksize, target_tablesize, num_tables)
    return sketch


def autoload(infile, count=True, graph=False, ksize=31, table_size=1e4,
             num_tables=4, num_bands=None, band=None):
    """
    Use file extension to conditionally load sketch into memory.

    If the file extension is one of the following, treat the file as a sketch
    that has been written to disk and load it with `kevlar.sketch.load`.
    Sketch attributes such as ksize, table size, and number of tables will
    be set automatically.

    - `.ct` or `.counttable`: `Counttable`
    - `.sct` or `.smallcounttable`: `SmallCounttable`
    - `.nt` or `.nodetable`: `Nodetable`
    - `.cg` or `.countgraph`: `Countgraph`
    - `.scg` or `.smallcountgraph`: `SmallCountgraph`
    - `.ng` or `.nodegraph`: `Nodegraph`

    Otherwise, a sketch will be created using the specified arguments and the
    input file will be treated as a FASTA/FASTQ file to be loaded with
    `.consume_seqfile` or `.consume_seqfile_banding`.
    """
    try:
        return load(infile)
    except KevlarSketchTypeError:
        sketch = allocate(ksize, table_size, num_tables, count=count,
                          graph=graph, smallcount=False)
        if num_bands:
            assert band >= 0 and band < num_bands
            sketch.consume_seqfile_banding(infile, num_bands, band)
        else:
            sketch.consume_seqfile(infile)
        return sketch
