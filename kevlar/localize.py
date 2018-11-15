#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
import kevlar
from kevlar.reference import bwa_align, ReferenceCutout
import pysam
import re
from subprocess import Popen, PIPE
import sys
from tempfile import TemporaryFile, NamedTemporaryFile


class KevlarRefrSeqNotFoundError(ValueError):
    """Raised if the reference sequence cannot be found."""
    pass


class Localizer(object):
    def __init__(self, seedsize, incl=None, excl=None):
        self._positions = defaultdict(list)
        self._seedsize = seedsize
        self.inclpattern = incl
        self.exclpattern = excl

    def __len__(self):
        return sum([
            len(self._positions[s]) for s in self._positions
            if not self.ignore_seqid(s)
        ])

    def ignore_seqid(self, seqid):
        """Determine whether a sequence should be discarded.

        This can be used to exclude calls from alternate/decoy sequences,
        organellar genomes, and so on.
        """
        include = True
        exclude = False
        if self.inclpattern:
            include = re.search(self.inclpattern, seqid) is not None
        if self.exclpattern:
            exclude = re.search(self.exclpattern, seqid) is not None
        return exclude or not include

    def add_seed_match(self, seqid, pos):
        """Store a seed match for later processing."""
        self._positions[seqid].append(pos)

    def get_cutouts(self, refrseqs=None, delta=0, clusterdist=1000):
        """Compute the reference target sequence(s) for the given seed matches.

        The cutout is defined as the span of all seed matches, plus `delta`
        nucleotides on either side. If two seed matches are separated by more
        than `clusterdist` nucleotides, they are split into different reference
        target sequences.
        """
        for seqid in sorted(self._positions):
            if self.ignore_seqid(seqid):
                continue
            matchpos = sorted(self._positions[seqid])
            assert len(matchpos) > 0
            if refrseqs and seqid not in refrseqs:
                raise KevlarRefrSeqNotFoundError(seqid)

            def new_cutout(cluster):
                startpos = max(cluster[0] - delta, 0)
                endpos = cluster[-1] + self._seedsize + delta
                subseq = None
                if refrseqs:
                    endpos = min(endpos, len(refrseqs[seqid]))
                    subseq = refrseqs[seqid][startpos:endpos]
                defline = '{:s}_{:d}-{:d}'.format(seqid, startpos, endpos)
                return ReferenceCutout(defline, subseq)

            if not clusterdist:
                yield new_cutout(matchpos)
                continue

            cluster = list()
            prevpos = None
            for nextpos in matchpos:
                if prevpos:
                    dist = nextpos - prevpos
                    if dist > clusterdist:
                        yield new_cutout(cluster)
                        cluster = list()
                cluster.append(nextpos)
                prevpos = nextpos
            yield new_cutout(cluster)


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
    print('[kevlar::localize]', message, file=logstream)
    seeds = set()
    for partition in partstream:
        contigs = list(partition)
        for contig in contigs:
            for seed in decompose_seeds(contig.sequence, seedsize):
                seeds.add(kevlar.revcommin(seed))
    n = 0
    for n, seed in enumerate(sorted(seeds)):
        print('>seed{}\n{}'.format(n, seed), file=seedstream)
    seedstream.flush()
    message = 'contigs decomposed into {} seeds'.format(n)
    print('[kevlar::localize]', message, file=logstream)


def get_seed_matches(seedfile, refrfile, seedsize=51, logstream=sys.stdout):
    """Determine the position of all seeds with a single system call to BWA."""
    print('[kevlar::localize] computing seed matches', file=logstream)
    bwa_cmd = 'bwa mem -k {k} -T {k} -a -c 5000 {idx} {seeds}'.format(
        k=seedsize, idx=refrfile, seeds=seedfile
    )
    bwa_args = bwa_cmd.split()
    seed_index = defaultdict(set)
    for seqid, start, end, seq in bwa_align(bwa_args, seqfilename=seedfile):
        minseq = kevlar.revcommin(seq)
        seed_index[minseq].add((seqid, start))
    message = 'found positions for {} seeds'.format(len(seed_index))
    print('[kevlar::localize]', message, file=logstream)
    return seed_index


def cutout(contigs, refrseqs, seed_matches, seedsize=51, delta=50,
           maxdiff=None, inclpattern=None, exclpattern=None, debug=False,
           logstream=sys.stdout):
    """Compute reference target sequences for a set of partitioned contigs.

    Partition by partition, decompose contigs into seeds, determine the genomic
    location of each seed, calculated the span of all seeds (plus some
    extension delta), and cut out that interval of the genome.
    """
    localizer = kevlar.localize.Localizer(
        seedsize, incl=inclpattern, excl=exclpattern
    )
    for contig in contigs:
        for seed in decompose_seeds(contig.sequence, seedsize):
            seed = kevlar.revcommin(seed)
            if seed not in seed_matches:
                if debug:  # pragma: no cover
                    message = 'WARNING: no position for seed {}'.format(seed)
                    print('[kevlar::localize]', message, file=logstream)
                continue
            for seqid, position in seed_matches[seed]:
                localizer.add_seed_match(seqid, position)
    if maxdiff is None:
        maxcontiglen = max([len(c.sequence) for c in contigs])
        maxdiff = maxcontiglen * 3

    cutter = localizer.get_cutouts(refrseqs=refrseqs, delta=delta,
                                   clusterdist=maxdiff)
    for gdna in cutter:
        yield gdna


def localize(partstream, refrfile, seedsize=51, delta=50, maxdiff=None,
             inclpattern=None, exclpattern=None, debug=False,
             logstream=sys.stdout):
    """Generator wrapper for the reference target cutout procedure."""
    partdata = list(partstream)
    partitions = [part for partid, part in partdata]
    partids = [partid for partid, part in partdata]
    message = 'loaded {} read partitions into memory'.format(len(partitions))
    print('[kevlar::localize]', message, file=logstream)
    kevlar.reference.autoindex(refrfile, logstream)

    with NamedTemporaryFile(mode='w', suffix='.contigs.fa') as seedfile:
        contigs_2_seeds(partitions, seedfile, seedsize=seedsize,
                        logstream=logstream)
        message = 'seeds written to "{}"'.format(seedfile.name)
        print('[kevlar::localize]', message, file=logstream)
        seed_matches = get_seed_matches(seedfile.name, refrfile,
                                        seedsize=seedsize, logstream=logstream)

    if len(seed_matches) == 0:
        message = 'WARNING: no reference matches'
        print('[kevlar::localize]', message, file=logstream)
        return

    message = 'loading reference sequences into memory'
    print('[kevlar::localize]', message, file=logstream)
    refrseqs = kevlar.seqio.parse_seq_dict(kevlar.open(refrfile, 'r'))

    message = 'computing the reference target sequence for each partition'
    print('[kevlar::localize]', message, file=logstream)
    ncutouts = 0
    progress_indicator = kevlar.ProgressIndicator(
        '[kevlar::localize]     computed targets for {counter} partitions',
        interval=100, breaks=[1000, 10000, 100000], logstream=logstream,
    )
    for partid, contiglist in partdata:
        progress_indicator.update()
        cutter = cutout(
            contiglist, refrseqs, seed_matches, seedsize=seedsize, delta=delta,
            maxdiff=maxdiff, inclpattern=inclpattern, exclpattern=exclpattern,
            debug=False, logstream=logstream
        )
        for gdna in cutter:
            ncutouts += 1
            yield partid, gdna
    if ncutouts == 0:
        message = 'WARNING: no reference matches'
        print('[kevlar::localize]', message, file=logstream)
        return


def main(args):
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(args.contigs, 'r'))
    if args.part_id:
        pstream = kevlar.parse_single_partition(contigstream, args.part_id)
    else:
        pstream = kevlar.parse_partitioned_reads(contigstream)
    outstream = kevlar.open(args.out, 'w')
    localizer = localize(
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
