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
from khmer import CyclicCounttable as Counttable
import kevlar


def load_sample_seqfile(seqfiles, ksize, memory, maxfpr=0.2,
                        mask=None, maskmaxabund=1, numbands=None, band=None,
                        outfile=None, numthreads=1, logfile=sys.stderr):
    """
    Compute k-mer abundances for the specified sequence input.

    Expected input is a list of one or more FASTA/FASTQ files corresponding
    to a single sample. A counttable is created and populated with abundances
    of all k-mers observed in the input. If `mask` is provided, only k-mers not
    present in the mask will be loaded.
    """
    sketch = Counttable(ksize, memory / 4, 4)
    numreads = 0
    for seqfile in seqfiles:
        print('[kevlar::count]      loading from', seqfile, file=logfile)
        parser = khmer.ReadParser(seqfile)
        threads = list()
        for _ in range(numthreads):
            if mask:
                if numbands:
                    thread = threading.Thread(
                        target=sketch.consume_seqfile_banding_with_mask,
                        args=(parser, numbands, band, mask, ),
                    )
                else:
                    thread = threading.Thread(
                        target=sketch.consume_seqfile_with_mask,
                        args=(parser, mask, ),
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
        if not outfile.endswith(('.ct', '.counttable')):
            outfile += '.counttable'
        sketch.save(outfile)
        message += ';\n    saved to "{:s}"'.format(outfile)
    print('[kevlar::count]    ', message, file=logfile)

    return sketch


def main(args):
    if (args.num_bands is None) is not (args.band is None):
        raise ValueError('Must specify --num-bands and --band together')
    myband = args.band - 1 if args.band else None

    timer = kevlar.Timer()
    timer.start()

    sketch = load_sample_seqfile(
        args.seqfile, args.ksize, args.memory, args.max_fpr,
        numbands=args.num_bands, band=myband, numthreads=args.threads,
        outfile=args.counttable, logfile=args.logfile
    )

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::count]', message, file=args.logfile)
