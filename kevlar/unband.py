#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 2018 Battelle National Biodefense Institute
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
import sys
from tempfile import TemporaryDirectory


def create_batch_files(numbatches, tempdir):
    batchfiles = list()
    for i in range(numbatches):
        tempfn = '{dir:s}/kevlar-unband-batch{batch:d}.augfastq.gz'.format(
            dir=tempdir, batch=i
        )
        tempfh = kevlar.open(tempfn, 'w')
        batchfiles.append(tempfh)
    return batchfiles


def write_records_to_batches(recordstream, batchfiles):
    numbatches = len(batchfiles)
    message = 'writing records to {:d} temp batch files'.format(numbatches)
    kevlar.plog('[kevlar::unband]', message)
    progress_indicator = kevlar.ProgressIndicator(
        '[kevlar::unband]     processed {counter} reads',
        interval=1e5, breaks=[1e6, 1e7],
    )
    for record in recordstream:
        progress_indicator.update()
        batch = hash(record.name) % numbatches
        fh = batchfiles[batch]
        kevlar.print_augmented_fastx(record, fh)


def resolve_batch(batchfile):
    reads = dict()
    filename = batchfile.name
    batchfile.close()
    batchfile = kevlar.open(filename, 'r')
    reader = kevlar.parse_augmented_fastx(batchfile)
    for read in reader:
        if read.name not in reads:
            reads[read.name] = read
            continue
        else:
            for ikmer in read.annotations:
                reads[read.name].annotations.append(ikmer)
    for readname in sorted(reads):
        read = reads[readname]
        read.annotations.sort(key=lambda k: k.offset)
        yield read
    batchfile.close()


def resolve_batches(batchfiles):
    numbatches = len(batchfiles)
    message = 'resolving duplicate reads in {:d} batches'.format(numbatches)
    kevlar.plog('[kevlar::unband]', message)
    for n, batchfile in enumerate(batchfiles):
        for read in resolve_batch(batchfile):
            yield read
        kevlar.plog('[kevlar::unband]     batch {:d} complete'.format(n))
    kevlar.plog('[kevlar::unband] Done!')


def unband(recordstream, numbatches=16):
    with TemporaryDirectory() as tempdir:
        batchfiles = create_batch_files(numbatches, tempdir)
        write_records_to_batches(recordstream, batchfiles)
        for read in resolve_batches(batchfiles):
            yield read


def main(args):
    outstream = kevlar.open(args.out, 'w')
    records = kevlar.seqio.afxstream(args.infile)
    for read in unband(records, args.n_batches):
        kevlar.print_augmented_fastx(read, outstream)
