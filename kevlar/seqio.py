#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
from itertools import combinations, product
from networkx import Graph, connected_components
from sys import stdout, stderr, exit
import re
import khmer
import screed
import kevlar


def parse_fasta(data):
    """
    Load sequences in Fasta format.

    This generator function yields a tuple containing a defline and a sequence
    for each record in the Fasta data. Stolen shamelessly from
    http://stackoverflow.com/a/7655072/459780.
    """
    name, seq = None, []
    for line in data:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))


def parse_seq_dict(data):
    """Load sequences from a Fasta file into a dictionary."""
    seqs = dict()
    for defline, sequence in parse_fasta(data):
        seqid = defline[1:].replace('\t', ' ').split(' ')[0]
        assert seqid not in seqs
        seqs[seqid] = sequence
    return seqs


def parse_augmented_fastx(instream):
    """
    Read augmented Fast[q|a] records into memory.

    The parsed records will have .name, .sequence, and .quality defined (unless
    it's augmented Fasta), as well as a list of interesting k-mers. See
    http://kevlar.readthedocs.io/en/latest/formats.html#augmented-sequences for
    more information.
    """
    record = None
    for line in instream:
        if line.startswith(('@', '>')):
            if record is not None:
                yield record
            readid = line[1:].strip()
            seq = next(instream).strip()
            if line.startswith('@'):
                _ = next(instream)
                qual = next(instream).strip()
                record = screed.Record(name=readid, sequence=seq, quality=qual,
                                       ikmers=list())
            else:
                record = screed.Record(name=readid, sequence=seq,
                                       ikmers=list())
        elif line.endswith('#\n'):
            offset = len(line) - len(line.lstrip())
            line = line.strip()[:-1]
            abundances = re.split('\s+', line)
            kmer = abundances.pop(0)
            abundances = [int(a) for a in abundances]
            ikmer = kevlar.KmerOfInterest(sequence=kmer, offset=offset,
                                          abund=abundances)
            record.ikmers.append(ikmer)
    if record is not None:
        yield record


def print_augmented_fastx(record, outstream=stdout):
    """Write augmented records out to an .augfast[q|a] file."""
    khmer.utils.write_record(record, outstream)
    for kmer in sorted(record.ikmers, key=lambda k: k.offset):
        abundstr = ' '.join([str(a) for a in kmer.abund])
        print(' ' * kmer.offset, kmer.sequence, ' ' * 10, abundstr, '#',
              sep='', file=outstream)


def afxstream(filelist):
    for infile in filelist:
        fh = kevlar.open(infile, 'r')
        for record in parse_augmented_fastx(fh):
            yield record


def load_reads_and_kmers(instream, logstream=None):
    """
    Load reads into lookup tables for convenient access.

    The first table is a dictionary of reads indexed by read name, and the
    second table is a dictionary of read sets indexed by an interesting k-mer.
    """
    reads = dict()
    kmers = defaultdict(set)
    for n, record in enumerate(kevlar.parse_augmented_fastx(instream), 1):
        if logstream and n % 10000 == 0:  # pragma: no cover
            print('[kevlar::seqio]    loaded {:d} reads'.format(n),
                  file=logstream)
        reads[record.name] = record
        for kmer in record.ikmers:
            kmerseq = kevlar.revcommin(kmer.sequence)
            kmers[kmerseq].add(record.name)
    return reads, kmers


class AnnotatedReadSet(object):
    """
    Data structure for de-duplicating reads and combining annotated k-mers.

    The `kevlar novel` command produces output in an "augmented Fastq" format,
    with "interesting" (potentially novel) k-mers annotated like so.

            @read1
            TTAACTCTAGATTAGGGGCGTGACTTAATAAGGTGTGGGCCTAAGCGTCT
            +
            IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                         AGGGGCGTGACTTAATAAG          8 0 0#
                           GGGCGTGACTTAATAAGGT          8 0 0#

    Each line ending with a # shows an interesting k-mer and its abundances
    in the case and control samples. If `kevlar novel` is run in "k-mer
    banding" mode (with the `--batch` flag), the same read will typically show
    up in multiple output files, with different "interesting" k-mers annotated
    in each. This data structure supports de-duplicating reads that appear in
    multiple augmented Fastq files and combining their annotated k-mers.
    """

    def __init__(self, ksize, abundmem):
        self._reads = dict()
        self._counts = khmer.Counttable(ksize, abundmem / 4, 4)
        self._readcounts = defaultdict(int)
        self._ikmercounts = defaultdict(int)

        self._masked = defaultdict(int)
        self._lowabund = defaultdict(int)
        self._valid = defaultdict(int)

        self._novalidkmers_count = 0

    def __len__(self):
        return len(self._reads)

    def __iter__(self):
        for readid in self._reads:
            record = self._reads[readid]
            if len(record.ikmers) == 0:
                continue
            yield record

    @property
    def masked(self):
        return len(self._masked), sum(self._masked.values())

    @property
    def lowabund(self):
        return len(self._lowabund), sum(self._lowabund.values())

    @property
    def valid(self):
        return len(self._valid), sum(self._valid.values())

    @property
    def discarded(self):
        return self._novalidkmers_count

    @property
    def distinct_reads(self):
        return len(self._readcounts)

    @property
    def read_instances(self):
        return sum(self._readcounts.values())

    @property
    def distinct_ikmers(self):
        return len(self._ikmercounts)

    @property
    def ikmer_instances(self):
        return sum(self._ikmercounts.values())

    def add(self, newrecord):
        if newrecord.name in self._reads:
            record = self._reads[newrecord.name]
            assert record.sequence == newrecord.sequence
            record.ikmers.extend(newrecord.ikmers)
        else:
            self._reads[newrecord.name] = newrecord
            self._counts.consume(newrecord.sequence)

        self._readcounts[newrecord.name] += 1
        for kmer in newrecord.ikmers:
            minkmer = kevlar.revcommin(kmer.sequence)
            self._ikmercounts[minkmer] += 1

    def validate(self, mask=None, minabund=5):
        for readid in self._reads:
            record = self._reads[readid]

            validated_kmers = list()
            for kmer in record.ikmers:
                kmerseq = kevlar.revcommin(kmer.sequence)
                if mask and mask.get(kmerseq) > 0:
                    self._masked[kmerseq] += 1
                elif self._counts.get(kmerseq) < minabund:
                    self._lowabund[kmerseq] += 1
                else:
                    kmer.abund[0] = self._counts.get(kmerseq)
                    validated_kmers.append(kmer)
                    self._valid[kmerseq] += 1
            record.ikmers = validated_kmers
            if len(validated_kmers) == 0:
                self._novalidkmers_count += 1
