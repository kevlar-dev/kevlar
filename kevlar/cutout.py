#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
from kevlar.reference import bwa_align
import pysam
import re
from subprocess import Popen, PIPE
import sys
from tempfile import TemporaryFile, NamedTemporaryFile


def decompose_seeds(seq, seedsize):
    """Break up a sequence into seeds.

    Seeds are essentially k-mers, but we use different terminology for here two
    reasons.

    1. These sequences are used specifically as anchors.
    2. We want to distinguish between k size (for k-mer counting) and seed size
       (for definining reference target sequences), which often have different
       optimal values.
    """
    for i in range(len(seq) - seedsize + 1):
        yield seq[i:i+seedsize]


def contigs_2_seeds(partstream, seedstream, seedsize=51, logstream=sys.stdout):
    """Convert a stream of partitioned contigs to seeds and write to a file."""
    message = 'decomposing contigs into seeds of length {}'.format(seedsize)
    print('[kevlar::cutout]', message, file=logstream)
    seeds = set()
    for partition in partstream:
        contigs = list(partition)
        for contig in contigs:
            for seed in decompose_seeds(contig.sequence, seedsize):
                seeds.add(kevlar.revcommin(seed))
    for n, seed in enumerate(sorted(seeds)):
        print('>seed{}\n{}'.format(n, seed), file=seedstream)
    seedstream.flush()


def get_seed_matches(seedfile, refrfile, seedsize=51, logstream=sys.stdout):
    """Determine the position of all seeds with a single system call to BWA."""
    print('[kevlar::cutout] computing seed matches', file=logstream)
    bwa_cmd = 'bwa mem -k {k} -T {k} -a -c 5000 {idx} {seeds}'.format(
        k=seedsize, idx=refrfile, seeds=seedfile
    )
    bwa_args = bwa_cmd.split()
    seed_index = dict()
    for seqid, start, end, seq in bwa_align(bwa_args, seqfilename=seedfile):
        seed_index[seq] = (seqid, start)
    return seed_index


def localize(contigs, refrseqs, seed_matches, seedsize=51, delta=50,
             maxdiff=None, inclpattern=None, exclpattern=None,
             debug=False, logstream=sys.stdout):
    """Compute reference target sequences for a set of partitioned contigs.

    Partition by partition, decompose contigs into seeds, determine the genomic
    location of each seed, calculated the span of all seeds (plus some
    extension delta), and cut out that interval of the genome.
    """
    localizer = kevlar.localize.Localizer(
        seedsize, delta=delta, incl=inclpattern, excl=exclpattern
    )
    for contig in contigs:
        for seed in decompose_seeds(contig.sequence, seedsize):
            seed = kevlar.revcommin(seed)
            if seed not in seed_matches:
                if debug:  # pragma: no cover
                    message = 'WARNING: no position for seed {}'.format(seed)
                    print('[kevlar::cutout]', message, file=logstream)
                continue
            seqid, position = seed_matches[seed]
            localizer.add_seed_match(seqid, position)
    if maxdiff is None:
        maxcontiglen = max([len(c.sequence) for c in contigs])
        maxdiff = maxcontiglen * 3

    for gdna in localizer.get_cutouts(refrseqs=refrseqs, clusterdist=maxdiff):
        yield gdna


def cutout(partstream, refrfile, seedsize=51, delta=50, maxdiff=None,
           inclpattern=None, exclpattern=None, debug=False,
           logstream=sys.stdout):
    """Generator wrapper for the reference target cutout procedure."""
    partdata = list(partstream)
    partitions = [part for partid, part in partdata]
    partids = [partid for partid, part in partdata]
    message = 'loaded {} read partitions into memory'.format(len(partitions))
    print('[kevlar::cutout]', message, file=logstream)
    kevlar.reference.autoindex(refrfile, logstream)

    with NamedTemporaryFile(mode='w', suffix='.contigs.fa') as seedfile:
        contigs_2_seeds(partitions, seedfile, seedsize=seedsize,
                        logstream=logstream)
        message = 'seeds written to "{}"'.format(seedfile.name)
        print('[kevlar::cutout]', message, file=logstream)
        seed_matches = get_seed_matches(seedfile.name, refrfile,
                                        seedsize=seedsize, logstream=logstream)

    if len(seed_matches) == 0:
        message = 'WARNING: no reference matches'
        print('[kevlar::cutout]', message, file=logstream)
        return

    message = 'loading reference sequences into memory'
    print('[kevlar::cutout]', message, file=logstream)
    refrseqs = kevlar.seqio.parse_seq_dict(kevlar.open(refrfile, 'r'))

    message = 'computing the reference target sequence for each partition'
    print('[kevlar::cutout]', message, file=logstream)
    for partid, contiglist in partdata:
        cutter = localize(
            contiglist, refrseqs, seed_matches, seedsize=seedsize, delta=delta,
            maxdiff=maxdiff, inclpattern=inclpattern, exclpattern=exclpattern,
            debug=True, logstream=logstream
        )
        for gdna in cutter:
            yield partid, gdna


def main(args):
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(args.contigs, 'r'))
    if args.part_id:
        pstream = kevlar.parse_single_partition(contigstream, args.part_id)
    else:
        pstream = kevlar.parse_partitioned_reads(contigstream)
    outstream = kevlar.open(args.out, 'w')
    localizer = cutout(
        pstream, args.refr, seedsize=args.seed_size, delta=args.delta,
        maxdiff=args.max_diff, inclpattern=args.include,
        exclpattern=args.exclude, logstream=args.logfile
    )
    for part, gdna in localizer:
        seqname = gdna.defline
        if part is not None:
            seqname += ' kvcc={}'.format(part)
        record = kevlar.sequence.Record(name=seqname, sequence=gdna.sequence)
        kevlar.sequence.write_record(record, outstream)
