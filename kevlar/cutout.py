#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
import pysam
import re
from subprocess import Popen, PIPE
import sys
from tempfile import TemporaryFile, NamedTemporaryFile


def decompose_seeds(seq, seedsize):
    for i in range(len(seq) - seedsize + 1):
        yield seq[i:i+seedsize]


def contigs_2_seeds(partstream, seedstream, seedsize=51, logstream=sys.stdout):
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
    print('[kevlar::cutout] computing seed matches', file=logstream)
    seed_index = dict()
    bwa_cmd = 'bwa mem -k {k} -T {k} -a -c 5000 {idx} {seeds}'.format(
        k=seedsize, idx=refrfile, seeds=seedfile
    )
    bwa_args = bwa_cmd.split()
    with TemporaryFile() as samfile:
        bwaproc = Popen(bwa_args, stdout=samfile, stderr=PIPE)
        stdout, stderr = bwaproc.communicate()
        if bwaproc.returncode != 0:
            print(stderr, file=sys.stderr)
            raise kevlar.reference.KevlarBWAError('problem running BWA')
        samfile.seek(0)
        sam = pysam.AlignmentFile(samfile, 'r')
        for record in sam:
            if record.is_unmapped:
                continue
            seqid = sam.get_reference_name(record.reference_id)
            seed_index[record.seq] = (seqid, record.reference_start)
    return seed_index


def localize(contigs, refrseqs, seed_matches, seedsize=51, delta=50,
             maxdiff=None, inclpattern=None, exclpattern=None,
             debug=False, logstream=sys.stdout):
    localizer = kevlar.localize.Localizer(
        seedsize, delta=delta, incl=inclpattern, excl=exclpattern
    )
    for contig in contigs:
        for seed in decompose_seeds(contig.sequence, seedsize):
            seed = kevlar.revcommin(seed)
            if seed not in seed_matches:
                if debug:
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
           inclpattern=None, exclpattern=None, logstream=sys.stdout):
    partitions = list(partstream)
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
    for contiglist in partitions:
        ccmatch = re.search(r'kvcc=(\d+)', contiglist[0].name)
        cc = ccmatch.group(1) if ccmatch else None
        cutter = localize(
            contiglist, refrseqs, seed_matches, seedsize=seedsize, delta=delta,
            maxdiff=maxdiff, inclpattern=inclpattern, exclpattern=exclpattern,
            logstream=logstream
        )
        for gdna in cutter:
            yield cc, gdna


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
