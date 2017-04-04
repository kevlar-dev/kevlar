#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict
from itertools import combinations, product
from networkx import Graph
from sys import stdout
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


def parse_augmented_fastq(instream):
    record = None
    for line in instream:
        if line.startswith('@'):
            if record is not None:
                yield record
            readid = line[1:].strip()
            seq = next(instream).strip()
            _ = next(instream)
            qual = next(instream).strip()
            record = screed.Record(name=readid, sequence=seq, quality=qual,
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


def print_augmented_fastq(record, outstream=stdout):
    khmer.utils.write_record(record, outstream)
    for kmer in sorted(record.ikmers, key=lambda k: k.offset):
        abundstr = ' '.join([str(a) for a in kmer.abund])
        print(' ' * kmer.offset, kmer.sequence, ' ' * 10, abundstr, '#',
              sep='', file=outstream)


class AnnotatedReadSet(object):
    """
    Data structure for de-duplicating reads and combining annotated k-mers.

    The `kevlar find` command produces output in an "augmented Fastq" format,
    with "interesting" (potentially novel) k-mers annotated like so.

            @read1
            TTAACTCTAGATTAGGGGCGTGACTTAATAAGGTGTGGGCCTAAGCGTCT
            +
            IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                         AGGGGCGTGACTTAATAAG          8 0 0#
                           GGGCGTGACTTAATAAGGT          8 0 0#

    Each line ending with a # shows an interesting k-mer and its abundances
    in the case and control samples. If `kevlar find` is run in "k-mer banding"
    mode (with the `--batch` flag), the same read will typically show up in
    multiple output files, with different "interesting" k-mers annotated in
    each. This data structure supports de-duplicating reads that appear in
    multiple augmented Fastq files and combining their annotated k-mers.
    """

    def __init__(self):
        self._reads = dict()

        self._masked = defaultdict(int)
        self._lowabund = defaultdict(int)
        self._valid = defaultdict(int)

        self._novalidkmers_count = 0
        self._contam_seqs = 0

    def __len__(self):
        return len(self._reads)

    def __iter__(self):
        for readid in self._reads:
            record = self._reads[readid]
            if len(record.ikmers) == 0:
                continue
            elif hasattr(record, 'contam') and record.contam is True:
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
    def contam(self):
        return self._contam_seqs

    def add(self, newrecord):
        if newrecord.name in self._reads:
            record = self._reads[newrecord.name]
            assert record.sequence == newrecord.sequence
            record.ikmers.extend(newrecord.ikmers)
        else:
            self._reads[newrecord.name] = newrecord

    def validate(self, counts, refr=None, contam=None, minabund=5):
        for readid in self._reads:
            record = self._reads[readid]
            if contam:
                medcount, _, _ = contam.get_median_count(record.sequence)
                if medcount > 0:
                    record.contam = True
                    self._contam_seqs += 1
                    continue

            validated_kmers = list()
            for kmer in record.ikmers:
                kmerseq = kevlar.revcommin(kmer.sequence)
                if refr and refr.get(kmerseq) > 0:
                    self._masked[kmerseq] += 1
                elif counts.get(kmerseq) < minabund:
                    self._lowabund[kmerseq] += 1
                else:
                    validated_kmers.append(kmer)
                    self._valid[kmerseq] += 1
            record.ikmers = validated_kmers
            if len(validated_kmers) == 0:
                self._novalidkmers_count += 1

    def group_reads_by_novel_kmers(self, upint=10000, logstream=None):
        reads_by_novel_kmer = defaultdict(set)
        for n, read_name in enumerate(self._reads):
            if logstream and n > 0 and n % upint == 0:
                print('    store reads by novel k-mers:', n, file=logstream)
            record = self._reads[read_name]
            for novel_kmer in record.ikmers:
                kmer_seq = novel_kmer.sequence
                kmer_seq_rc = kevlar.revcom(novel_kmer.sequence)
                reads_by_novel_kmer[kmer_seq].add(record.name)
                reads_by_novel_kmer[kmer_seq_rc].add(record.name)

        read_graph = Graph()
        for n, read_name in enumerate(self._reads):
            if logstream and n > 0 and n % upint == 0:
                print('    build shared novel k-mer graph:', n, file=logstream)
            record = self._reads[read_name]
            for kmer in record.ikmers:
                for other_record_name in reads_by_novel_kmer[kmer.sequence]:
                    read_graph.add_edge(read_name, other_record_name)

        for cc in read_graph.connected_components():
            print('DEBUG', cc, file=sys.stderr)
