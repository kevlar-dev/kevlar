#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import khmer
from khmer import khmer_args
import kevlar
from kevlar.count import load_sample_seqfile
from kevlar.sketch import KevlarUnsuitableFPRError
from kevlar.sketch import sketch_loader_by_filename_extension as sketch_ext


def load_mask(maskfiles, ksize, memory, maxfpr=0.001, savefile=None,
              numthreads=1):
    if maskfiles is None:
        return None
    if len(maskfiles) == 1 and maskfiles[0].endswith(tuple(sketch_ext)):
        mask = kevlar.sketch.load(maskfiles[0])
        kevlar.plog('[kevlar::filter] Mask "{}" loaded'.format(maskfiles[0]))
    else:
        mask = load_sample_seqfile(
            maskfiles, ksize, memory, maxfpr=maxfpr, count=False,
            outfile=savefile, numthreads=numthreads
        )
    return mask


def summarize_readset(readset):
    fpr = kevlar.sketch.estimate_fpr(readset._counts)
    message = '    {:d} instances'.format(readset.read_instances)
    message += ' of {:d} reads consumed,\n'.format(readset.distinct_reads)
    message += '    annotated with'
    message += ' {:d} instances '.format(readset.ikmer_instances)
    message += 'of {:d} distinct'.format(readset.distinct_ikmers)
    message += ' "interesting" k-mers;\n'
    message += '    estimated false positive rate is {:1.3f}'.format(fpr)
    kevlar.plog(message)
    return fpr


def summarize_validate(readset, n):
    int_distinct = readset.masked[0] + readset.abundfilt[0] + readset.valid[0]
    int_instances = readset.masked[1] + readset.abundfilt[1] + readset.valid[1]

    message = '    processed {:d} instances'.format(int_instances)
    message += ' of {:d} distinct "interesting" k-mers'.format(int_distinct)
    message += ' in {:d} reads'.format(len(readset))
    message += '\n        '
    message += '{:d} instances'.format(readset.masked[1])
    message += ' of {:d} distinct k-mers'.format(readset.masked[0])
    message += ' masked as contaminant or included in reference genome'
    message += '\n        '
    message += '{:d} instances'.format(readset.abundfilt[1])
    message += ' of {:d} distinct k-mers'.format(readset.abundfilt[0])
    message += ' discarded due to abundance filtering'
    message += '\n        '
    message += '{:d} instances'.format(readset.valid[1])
    message += ' of {:d} distinct k-mers'.format(readset.valid[0])
    message += ' validated as novel'
    message += '\n        '
    message += '{:d} reads'.format(readset.discarded)
    message += ' with no surviving valid k-mers ignored'
    message += '\n        '
    message += '{:d} reads written to output'.format(n)
    kevlar.plog(message)


def filter(readstream, mask=None, casemin=6, ctrlmax=1, ksize=31, memory=1e6,
           maxfpr=0.001):
    timer = kevlar.Timer()
    timer.start('recalc')
    kevlar.plog(
        '[kevlar::filter] Loading input; recalculate k-mer abundances,',
        'de-duplicate reads and merge k-mers',
    )
    readset = kevlar.seqio.AnnotatedReadSet(ksize, memory, mask=mask)
    for record in readstream:
        readset.add(record)
    fpr = summarize_readset(readset)
    if fpr > maxfpr:  # pragma: no cover
        raise KevlarUnsuitableFPRError('FPR too high, bailing out!!!')
    elapsed = timer.stop('recalc')
    kevlar.plog('[kevlar::filter] Input loaded in {:.2f} sec'.format(elapsed))

    timer.start('validate')
    kevlar.plog('[kevlar::filter] Validate k-mers and print reads')
    readset.validate(casemin=casemin, ctrlmax=ctrlmax)
    n = 0
    for n, record in enumerate(readset, 1):
        yield record
    summarize_validate(readset, n)
    elapsed = timer.stop('validate')
    kevlar.plog(
        '[kevlar::filter] k-mers validated and reads printed',
        'in {:.2f} sec'.format(elapsed),
    )


def main(args):
    timer = kevlar.Timer()
    timer.start()

    mask = None
    if args.mask:
        message = 'Loading mask (sequences to ignore/discard)'
        kevlar.plog('[kevlar::filter]', message)
        mask = load_mask(
            args.mask, args.ksize, args.mask_memory, maxfpr=args.mask_max_fpr,
            savefile=args.save_mask, numthreads=args.threads,
        )

    readstream = kevlar.seqio.afxstream(args.augfastq)
    outstream = kevlar.open(args.out, 'w')
    filterstream = filter(
        readstream, mask, casemin=args.case_min, ctrlmax=args.ctrl_max,
        ksize=args.ksize, memory=args.abund_memory, maxfpr=args.abund_max_fpr,
    )
    for record in filterstream:
        kevlar.print_augmented_fastx(record, outstream)

    total = timer.stop()
    message = 'Total time: {:.2f} seconds'.format(total)
    kevlar.plog('[kevlar::filter]', message)
