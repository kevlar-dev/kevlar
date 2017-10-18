#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import re
import sys

import khmer
from khmer import khmer_args
import kevlar


def load_mask(maskfiles, ksize, memory, maxfpr=0.001, savefile=None,
              logfile=sys.stderr):
    """Load reference genome and/or contaminant database from a file."""
    if len(maskfiles) == 1 and maskfiles[0].endswith(('.nt', '.nodetable')):
        mask = kevlar.sketch.load(maskfiles[0])
        message = '    nodetable loaded'
    else:
        buckets = memory * khmer._buckets_per_byte['nodegraph'] / 4
        mask = khmer.Nodetable(ksize, buckets, 4)
        nr, nk = 0, 0
        for maskfile in maskfiles:
            numreads, numkmers = mask.consume_seqfile(maskfile)
            nr += numreads
            nk += numkmers
        message = '    {:d} sequences and {:d} k-mers consumed'.format(nr, nk)
    fpr = kevlar.sketch.estimate_fpr(mask)
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    print(message, file=logfile)
    if fpr > maxfpr:
        print('[kevlar::filter] FPR too high, bailing out', file=logfile)
        sys.exit(1)
    if savefile:
        mask.save(savefile)
        message = '    nodetable saved to "{:s}"'.format(savefile)
        print(message, file=logfile)
    return mask


def load_input(filelist, ksize, memory, maxfpr=0.001, logfile=sys.stderr):
    """
    Load input data.

    The input data is loaded into two data structures. First, the read
    sequences are loaded into a countgraph to recompute k-mer abundances with
    (effectively) exact precision. Second, the reads and their corresponding
    "interesting" k-mers are loaded into an AnnotatedReadSet to de-duplicate
    reads and group k-mers by read.
    """
    countgraph = khmer.Countgraph(ksize, memory / 4, 4)
    read_inst_consumed = 0
    int_kmer_instances = 0
    int_kmers_parsed = set()
    readset = kevlar.seqio.AnnotatedReadSet()
    for filename in filelist:
        print('    -', filename, file=logfile)
        with kevlar.open(filename, 'r') as infile:
            for record in kevlar.parse_augmented_fastx(infile):
                if record.name not in readset._reads:
                    countgraph.consume(record.sequence)
                readset.add(record)
                read_inst_consumed += 1
                for kmer in record.ikmers:
                    int_kmer_instances += 1
                    minkmer = kevlar.revcommin(kmer.sequence)
                    int_kmers_parsed.add(minkmer)
    n_kmers_distinct = len(int_kmers_parsed)

    fpr = kevlar.sketch.estimate_fpr(countgraph)
    message = '    {:d} instances'.format(read_inst_consumed)
    message += ' of {:d} reads consumed'.format(len(readset))
    message += ', annotated with {:d} instances '.format(int_kmer_instances)
    message += 'of {:d} distinct "interesting" k-mers'.format(n_kmers_distinct)
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    print(message, file=logfile)
    if fpr > maxfpr:
        print('[kevlar::filter] FPR too high, bailing out', file=logfile)
        sys.exit(1)
    return readset, countgraph


def validate_and_print(readset, countgraph, mask=None, minabund=5,
                       outfile=sys.stdout, logfile=sys.stderr):
    readset.validate(countgraph, mask=mask, minabund=minabund)
    n = 0  # Get an unbound var error later (printing report) without this?!?!
    for n, record in enumerate(readset):
        kevlar.print_augmented_fastx(record, outfile)

    int_distinct = readset.masked[0] + readset.lowabund[0] + readset.valid[0]
    int_instances = readset.masked[1] + readset.lowabund[1] + readset.valid[1]

    message = '    processed {:d} instances'.format(int_instances)
    message += ' of {:d} distinct "interesting" k-mers'.format(int_distinct)
    message += ' in {:d} reads'.format(len(readset))
    message += '\n        '
    message += '{:d} instances'.format(readset.masked[1])
    message += ' of {:d} distinct k-mers'.format(readset.masked[0])
    message += ' masked by the reference genome'
    message += '\n        '
    message += '{:d} instances'.format(readset.lowabund[1])
    message += ' of {:d} distinct k-mers'.format(readset.lowabund[0])
    message += ' discarded due to low abundance'
    message += '\n        '
    message += '{:d} instances'.format(readset.valid[1])
    message += ' of {:d} distinct k-mers'.format(readset.valid[0])
    message += ' validated as novel'
    message += '\n        '
    message += '{:d} reads'.format(readset.discarded)
    message += ' with no surviving valid k-mers ignored'
    message += '\n        '
    message += '{:d} reads written to output'.format(n + 1)
    print(message, file=logfile)


def main(args):
    timer = kevlar.Timer()
    timer.start()

    mask = None
    if args.mask:
        timer.start('loadmask')
        print('[kevlar::filter] Loading mask from', args.mask,
              file=args.logfile)
        mask = load_mask(args.mask, args.ksize, args.mask_memory,
                         maxfpr=args.mask_max_fpr, savefile=args.save_mask,
                         logfile=args.logfile)
        elapsed = timer.stop('loadmask')
        print('[kevlar::filter]', 'Mask loaded in {:.2f} sec'.format(elapsed),
              file=args.logfile)

    timer.start('recalc')
    print('[kevlar::filter] Loading input; recalculate k-mer abundances,',
          'de-duplicate reads and merge k-mers',
          file=args.logfile)
    readset, countgraph = load_input(args.augfastq, args.ksize,
                                     args.abund_memory, args.abund_max_fpr,
                                     args.logfile)
    elapsed = timer.stop('recalc')
    print('[kevlar::filter] Input loaded in {:.2f} sec'.format(elapsed),
          file=args.logfile)

    timer.start('validate')
    print('[kevlar::filter] Validate k-mers and print reads',
          file=args.logfile)
    outstream = kevlar.open(args.out, 'w')
    validate_and_print(readset, countgraph, mask, args.min_abund, outstream,
                       args.logfile)
    elapsed = timer.stop('validate')
    print('[kevlar::filter] k-mers validated and reads printed',
          'in {:.2f} sec'.format(elapsed), file=args.logfile)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    print('[kevlar::filter]', message, file=args.logfile)
