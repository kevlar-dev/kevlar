#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import khmer
import sys


sketch_loader_by_filename_extension = {
    '.nt':  khmer.Nodetable.load,
    '.ng':  khmer.Nodegraph.load,
    '.ct':  khmer.Counttable.load,
    '.cg':  khmer.Countgraph.load,
    '.sct': khmer.SmallCounttable.load,
    '.scg': khmer.SmallCountgraph.load,
    '.nodetable':       khmer.Nodetable.load,
    '.nodegraph':       khmer.Nodegraph.load,
    '.counttable':      khmer.Counttable.load,
    '.countgraph':      khmer.Countgraph.load,
    '.smallcounttable': khmer.SmallCounttable.load,
    '.smallcountgraph': khmer.SmallCountgraph.load,
}

# count(?)->graph(?)->small(?)
sketch_extensions_by_trait = {
    True: {
        True: {
            True: ('.scg', '.smallcountgraph'),
            False: ('.cg', '.countgraph'),
        },
        False: {
            True: ('.sct', '.smallcounttable'),
            False: ('.ct', '.counttable'),
        },
    },
    False: {
        True: {
            True: ('.ng', '.nodegraph'),
            False: ('.ng', '.nodegraph'),
        },
        False: {
            True: ('.nt', '.nodetable'),
            False: ('.nt', '.nodetable'),
        },
    },
}


class KevlarSketchTypeError(ValueError):
    pass


class KevlarUnsuitableFPRError(SystemExit):
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
    extensions = tuple(sketch_loader_by_filename_extension)
    if not filename.endswith(extensions):
        message = 'unable to determine sketch type from filename ' + filename
        raise KevlarSketchTypeError(message)
    ext = '.' + filename.split('.')[-1]
    loadfunc = sketch_loader_by_filename_extension[ext]
    return loadfunc(filename)


def get_extension(count=False, graph=False, smallcount=False):
    return sketch_extensions_by_trait[count][graph][smallcount]


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


def load_sketchfiles(sketchfiles, maxfpr=0.2, logfile=sys.stderr):
    """Load samples from pre-computed k-mer abundances."""
    sketches = list()
    for sketchfile in sketchfiles:
        message = 'loading sketchfile "{}"...'.format(sketchfile)
        print('[kevlar::sketch]    ', message, end='', file=logfile)
        sketch = autoload(sketchfile)
        fpr = estimate_fpr(sketch)
        message = 'done! estimated false positive rate is {:1.3f}'.format(fpr)
        if fpr > maxfpr:
            message += ' (FPR too high, bailing out!!!)'
            raise KevlarUnsuitableFPRError(message)
        print(message, file=logfile)
        sketches.append(sketch)
    return sketches
