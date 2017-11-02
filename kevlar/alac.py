#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import kevlar
from kevlar.assemble import assemble_greedy, assemble_fml_asm
from kevlar.localize import localize
from kevlar.call import call


def alac(readstream, refrfile, ksize=31, delta=25, maxdiff=10000, match=1,
         mismatch=2, gapopen=5, gapextend=0, greedy=False,
         logstream=sys.stderr):
    assembler = assemble_greedy if greedy else assemble_fml_asm
    cntgs = [c for c in assembler(readstream, logstream=logstream)]
    targets = [t for t in localize(cntgs, refrfile, ksize=ksize, delta=delta)]
    caller = call(targets, cntgs, match, mismatch, gapopen, gapextend, ksize)
    for varcall in caller:
        yield varcall


def main(args):
    readstream = kevlar.parse_augmented_fastx(kevlar.open(args.infile, 'r'))
    outstream = kevlar.open(args.out, 'w')
    workflow = alac(
        readstream, args.refr, ksize=args.ksize, delta=args.delta,
        maxdiff=args.max_diff, match=args.match, mismatch=args.mismatch,
        gapopen=args.open, gapextend=args.extend, greedy=args.greedy,
        logstream=args.logfile
    )

    for varcall in workflow:
        print(varcall.vcf, file=outstream)
