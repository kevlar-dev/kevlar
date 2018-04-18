#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from itertools import chain
from queue import Queue
import re
import sys
from threading import Thread
from time import sleep

import kevlar
from kevlar.assemble import assemble_greedy, assemble_fml_asm
from kevlar.localize import localize
from kevlar.call import call


def populate_queue(queue, idx, pstream, bigpart=10000, logstream=sys.stderr):
    for partition in pstream:
        reads = list(partition)
        if len(reads) > bigpart:
            message = '[kevlar::alac::populate_queue] WARNING: '
            message += 'skipping partition with {:d} reads'.format(len(reads))
            print(message, file=sys.stderr)
            continue
        message = '[kevlar::alac::populate_queue (thread={})] '.format(idx)
        message += 'adding partition with {} reads to queue'.format(len(reads))
        print(message, file=sys.stderr)
        queue.put(reads)


def make_call_from_reads(queue, idx, calls, refrfile, ksize=31, delta=50,
                         seedsize=31, maxdiff=None, match=1, mismatch=2,
                         gapopen=5, gapextend=0, greedy=False, fallback=False,
                         min_ikmers=None, refrseqs=None, logstream=sys.stderr):
        while True:
            if queue.empty():
                sleep(1)
                continue
            reads = queue.get()
            ccmatch = re.search(r'kvcc=(\d+)', reads[0].name)
            cc = ccmatch.group(1) if ccmatch else None
            message = '[kevlar::alac::make_call_from_reads'
            message += ' (thread={:d})]'.format(idx)
            message += ' grabbed partition={} from queue,'.format(cc)
            message += ' queue size now {:d}'.format(queue.qsize())
            print(message, file=sys.stderr)

            # Assemble partitioned reads into contig(s)
            assmblr = assemble_greedy if greedy else assemble_fml_asm
            contigs = list(assmblr(reads, logstream=logstream))
            if len(contigs) == 0 and assmblr == assemble_fml_asm and fallback:
                message = 'WARNING: no contig assembled by fermi-lite'
                if cc:
                    message += ' for partition={:s}'.format(cc)
                message += '; attempting again with homegrown greedy assembler'
                print('[kevlar::alac]', message, file=logstream)
                contigs = list(assemble_greedy(reads, logstream=logstream))
            if min_ikmers is not None:
                # Apply min ikmer filter if it's set
                contigs = [c for c in contigs if len(c.ikmers) >= min_ikmers]
            if len(contigs) == 0:
                queue.task_done()
                continue

            # Identify the genomic region(s) associated with each contig
            localizer = localize(
                contigs, refrfile, seedsize, delta=delta, maxdiff=maxdiff,
                refrseqs=refrseqs, logstream=logstream
            )
            targets = list(localizer)
            if len(targets) == 0:
                queue.task_done()
                continue

            # Align contigs to genomic targets to make variant calls
            caller = call(targets, contigs, match, mismatch, gapopen,
                          gapextend, ksize, refrfile)
            for varcall in caller:
                if cc:
                    varcall.annotate('PART', cc)
                calls.append(varcall)
            queue.task_done()


def alac(pstream, refrfile, threads=1, ksize=31, bigpart=10000, delta=50,
         seedsize=31, maxdiff=None, match=1, mismatch=2, gapopen=5,
         gapextend=0, greedy=False, fallback=False, min_ikmers=None,
         logstream=sys.stderr):
    part_queue = Queue(maxsize=max(32, 12 * threads))

    refrstream = kevlar.open(refrfile, 'r')
    refrseqs = kevlar.seqio.parse_seq_dict(refrstream)

    call_lists = list()
    for idx in range(threads + 1):
        if idx == 0:
            worker = Thread(
                target=populate_queue,
                args=(part_queue, idx, pstream, bigpart, logstream,),
            )
            worker.start()
            continue
        thread_calls = list()
        call_lists.append(thread_calls)
        worker = Thread(
            target=make_call_from_reads,
            args=(
                part_queue, idx, thread_calls, refrfile, ksize, delta,
                seedsize, maxdiff, match, mismatch, gapopen, gapextend, greedy,
                fallback, min_ikmers, refrseqs, logstream,
            )
        )
        worker.setDaemon(True)
        worker.start()

    part_queue.join()
    allcalls = sorted(chain(*call_lists), key=lambda c: (c.seqid, c.position))
    for call in allcalls:
        yield call


def main(args):
    readstream = kevlar.parse_augmented_fastx(kevlar.open(args.infile, 'r'))
    if args.part_id:
        pstream = kevlar.parse_single_partition(readstream, args.part_id)
    else:
        pstream = kevlar.parse_partitioned_reads(readstream)
    outstream = kevlar.open(args.out, 'w')
    workflow = alac(
        pstream, args.refr, threads=args.threads, ksize=args.ksize,
        bigpart=args.bigpart, delta=args.delta, seedsize=args.seed_size,
        maxdiff=args.max_diff, match=args.match, mismatch=args.mismatch,
        gapopen=args.open, gapextend=args.extend, greedy=args.greedy,
        fallback=args.fallback, min_ikmers=args.min_ikmers,
        logstream=args.logfile
    )

    for varcall in workflow:
        print(varcall.vcf, file=outstream)
