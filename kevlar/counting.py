#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import re
import sys

import khmer
import kevlar


def table_filename(seqfile, band=None):
    """
    asdf
    """
    suffixes = ['.fa', '.fasta', '.fna', '.fq', '.fastq']
    gzipsuffixes = [s + '.gz' for s in suffixes]
    if seqfile.endswith(tuple(suffixes)):
        prefix = '.'.join(seqfile.split('.')[:-1])
    elif sample.endswith(tuple(gzipsuffixes)):
        prefix = '.'.join(seqfile.split('.')[:-2])
    else:
        prefix = seqfile
    if band:
        prefix += '.band{:d}'.format(band)
    return prefix + '.counttable'


def load_sample_seqfile(seqfiles, ksize, memory, maxfpr=0.2, masks=None,
                        maskmaxabund=1, numbands=None, band=None,
                        skipsave=False, logfile=sys.stderr):
    """
    asdf
    """
    message = 'loading sample from ' + ','.join(seqfiles)
    print('[kevlar::counting]    ', message, file=logfile)

    sketch = khmer.Counttable(ksize, memory / 4, 4)
    n, nkmers = 0, 0
    for n, read in kevlar.multi_file_iter_khmer(seqfiles, 1):
        cleansubseqs = kevlar.clean_subseqs(read.sequence)
        for kmer in sketch.get_kmers(cleansubseqs):
            if numbands:
                khash = sketch.hash(kmer)
                if sketch & (numbands - 1) != band - 1:
                    continue
            if masks:
                for mask in masks:
                    if mask.get(kmer) > maskmaxabund:
                        break
                else:
                    sketch.add(kmer)
                    nkmers += 1

    message = 'done loading reads'
    if numbands:
        message += ' (band {:d}/{:d})'.format(band, numbands)
    fpr = kevlar.calc_fpr(ct)
    message += '; {:d} reads processed'.format(n)
    message += ', {:d} k-mers stored'.format(nkmers)
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    if fpr > max_fpr:
        message += ' (FPR too high, bailing out!!!)'
        raise SystemExit(message)
    else:
        if not skipsave:
            savename = table_filename(sample, band)
            sketch.save(savename)
            message += '; saved to "{:s}"'.format(savename)
        print('[kevlar::counting]    ', message, file=logfile)

    return sketch


def load_samples_with_dilution(seqfilelists, ksize, memory, memfraction=0.1,
                               maxfpr=0.2, maxabund=1, numbands=None,
                               band=None, skipsave=False, logfile=sys.stderr):
    """
    asdf
    """
    numsamples = len(seqfilelists)
    message = 'computing k-mer abundances for {:d} samples'.format(numsamples)
    print('[kevlar::counting]    ', message, file=logfile)

    sketches = list()
    for seqfiles in seqfilelists:
        if len(sketches) == 0:
            masks = None
            sketchmem = memory
        else:
            masks = sketches
            sketchmem = memory * memfraction
        sketch = load_sample_seqfile(
            seqfiles, ksize, sketchmem, maxfpr=maxfpr, masks=masks,
            maskmaxabund=maxabund, numbands=numbands, band=band,
            skipsave=skipsave, logfile=logfile
        )
        sketches.append(sketch)
    return sketches


def load_samples_sketchfiles(sketchfiles, maxfpr=0.2, logfile=sys.stderr):
    """
    asdf
    """
    sketches = list()
    for sketchfile in sketchfiles:
        message = 'loading sketchfile ' + sketchfile,
        print('[kevlar::counting]    ', message, end='', file=logfile)
        sketch = kevlar.sketch_autoload(sketchfile)
        fpr = kevlar.calc_fpr(sketch)
        message = 'done! estimated false positive rate is {:1.3f}'.format(fpr)
        if fpr > maxfpr:
            message += ' (FPR too high, bailing out!!!)'
            raise SystemExit(message)
        else:
            print(message, file=logfile)
        sketches.append(sketch)
    return sketches
