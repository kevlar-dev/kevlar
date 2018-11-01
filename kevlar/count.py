#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
import threading
import sys
import khmer
import kevlar
from kevlar.sketch import allocate, get_extension


def load_sample_seqfile(seqfiles, ksize, memory, maxfpr=0.2, count=True,
                        smallcount=False, mask=None, maskmaxabund=0,
                        consume_masked=False, numbands=None, band=None,
                        outfile=None, numthreads=1, logfile=sys.stderr):
    """Compute k-mer abundances for the specified sequence input.

    Expected input is a list of one or more FASTA/FASTQ files corresponding
    to a single sample. A sketch is created and populated with abundances
    of all k-mers observed in the input. If `mask` is provided, only k-mers not
    present in the mask will be loaded.
    """
    numtables = 4
    sketchtype = 'nodegraph'
    if count:
        sketchtype = 'smallcountgraph' if smallcount else 'countgraph'
    tablesize = memory / numtables * khmer._buckets_per_byte[sketchtype]
    sketch = allocate(ksize, tablesize, num_tables=numtables, count=count,
                      smallcount=smallcount)
    numreads = 0
    for seqfile in seqfiles:
        print('[kevlar::count]      loading from', seqfile, file=logfile)
        parser = khmer.ReadParser(seqfile)
        threads = list()
        for _ in range(numthreads):
            if mask:
                threshold = 1 if consume_masked else maskmaxabund
                kwargs = {
                    'consume_masked': consume_masked,
                    'threshold': threshold
                }
                if numbands:
                    thread = threading.Thread(
                        target=sketch.consume_seqfile_banding_with_mask,
                        args=(parser, numbands, band, mask, ),
                        kwargs=kwargs,
                    )
                else:
                    thread = threading.Thread(
                        target=sketch.consume_seqfile_with_mask,
                        args=(parser, mask, ),
                        kwargs=kwargs,
                    )
            else:
                if numbands:
                    thread = threading.Thread(
                        target=sketch.consume_seqfile_banding,
                        args=(parser, numbands, band, ),
                    )
                else:
                    thread = threading.Thread(
                        target=sketch.consume_seqfile,
                        args=(parser, ),
                    )
            threads.append(thread)
            thread.start()

        for thread in threads:
            thread.join()
        numreads += parser.num_reads

    message = 'done loading reads'
    if numbands:
        message += ' (band {:d}/{:d})'.format(band+1, numbands)
    fpr = kevlar.sketch.estimate_fpr(sketch)
    message += ';\n    {:d} reads processed'.format(numreads)
    message += ', {:d} distinct k-mers stored'.format(sketch.n_unique_kmers())
    message += ';\n    estimated false positive rate is {:1.3f}'.format(fpr)
    if fpr > maxfpr:
        message += ' (FPR too high, bailing out!!!)'
        message = '[kevlar::count]     ' + message
        raise kevlar.sketch.KevlarUnsuitableFPRError(message)

    if outfile:
        extensions = get_extension(count=count, smallcount=smallcount)
        if not outfile.endswith(extensions):
            outfile += extensions[1]
        sketch.save(outfile)
        message += ';\n    saved to "{:s}"'.format(outfile)
    print('[kevlar::count]    ', message, file=logfile)

    return sketch


def main(args):
    if (args.num_bands is None) is not (args.band is None):
        raise ValueError('Must specify --num-bands and --band together')
    myband = args.band - 1 if args.band else None
    if args.mask:
        args.mask = kevlar.sketch.load(args.mask)

    timer = kevlar.Timer()
    timer.start()

    docount = args.counter_size > 1
    dosmallcount = args.counter_size == 4
    sketch = load_sample_seqfile(
        args.seqfile, args.ksize, args.memory, args.max_fpr, count=docount,
        smallcount=dosmallcount, mask=args.mask,
        consume_masked=args.count_masked, numbands=args.num_bands, band=myband,
        numthreads=args.threads, outfile=args.counttable, logfile=args.logfile
    )

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::count]', message, file=args.logfile)
