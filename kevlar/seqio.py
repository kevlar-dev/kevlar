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
from kevlar.sequence import Record, KmerOfInterest
from kevlar.sequence import write_record, parse_augmented_fastx


class KevlarPartitionLabelError(ValueError):
    pass


def parse_fasta(data):
    """Load sequences in Fasta format.

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
    if name:  # pragma: no cover
        yield (name, ''.join(seq))


def parse_seq_dict(data):
    """Load sequences from a Fasta file into a dictionary."""
    seqs = dict()
    for defline, sequence in parse_fasta(data):
        seqid = defline[1:].replace('\t', ' ').split(' ')[0]
        assert seqid not in seqs, seqid
        seqs[seqid] = sequence
    return seqs


def afxstream(filelist):
    for infile in filelist:
        fh = kevlar.open(infile, 'r')
        for record in parse_augmented_fastx(fh):
            yield record


def partition_id(readname):
    partmatch = re.search(r'kvcc=(\d+)', readname)
    if not partmatch:
        return None
    return partmatch.group(1)


def parse_partitioned_reads(readstream):
    current_part = None
    reads = list()
    for read in readstream:
        name = read.name if hasattr(read, 'name') else read.defline
        part = partition_id(name)
        if part is None:
            reads.append(read)
            current_part = False
            continue

        if current_part is False:
            message = 'reads with and without partition labels (kvcc=#)'
            raise KevlarPartitionLabelError(message)

        if part != current_part:
            if current_part:
                yield current_part, reads
                reads = list()
            current_part = part
        reads.append(read)

    if current_part is False:
        current_part = None
    yield current_part, reads


def parse_single_partition(readstream, partid):
    """
    Retrieve a single partition (by label) from a stream of partitioned reads.
    """
    for pid, partition in parse_partitioned_reads(readstream):
        if pid == partid:
            yield pid, partition


class AnnotatedReadSet(object):
    """Data structure for de-duplicating reads and combining annotated k-mers.

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

    def __init__(self, ksize, abundmem, mask=None):
        self._reads = dict()
        self._counts = khmer.Counttable(ksize, abundmem / 4, 4)
        self._mask = mask
        self._readcounts = defaultdict(int)
        self._ikmercounts = defaultdict(int)

        self._masked = defaultdict(int)
        self._abundfilt = defaultdict(int)
        self._valid = defaultdict(int)

        self._novalidkmers_count = 0

    def __len__(self):
        return len(self._reads)

    def __iter__(self):
        for readid in self._reads:
            record = self._reads[readid]
            if len(record.annotations) == 0:
                continue
            yield record

    @property
    def masked(self):
        return len(self._masked), sum(self._masked.values())

    @property
    def abundfilt(self):
        return len(self._abundfilt), sum(self._abundfilt.values())

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
        self._readcounts[newrecord.name] += 1
        if newrecord.name in self._reads:
            record = self._reads[newrecord.name]
            assert record.sequence == newrecord.sequence
            record.annotations.extend(newrecord.annotations)
        else:
            self._reads[newrecord.name] = newrecord

        for kmer in newrecord.annotations:
            kmerhash = self._counts.hash(newrecord.ikmerseq(kmer))
            self._ikmercounts[kmerhash] += 1
            if self._mask and self._mask.get(kmerhash) > 0:
                self._masked[kmerhash] += 1
            else:
                self._counts.add(kmerhash)

    def validate(self, casemin=6, ctrlmax=1):
        for readid in self._reads:
            record = self._reads[readid]
            validated_kmers = list()
            for kmer in record.annotations:
                kmerhash = self._counts.hash(record.ikmerseq(kmer))
                kmercount = self._counts.get(kmerhash)
                if kmercount < casemin:
                    self._abundfilt[kmerhash] += 1
                elif sum([1 for a in kmer.abund[1:] if a > ctrlmax]):
                    self._abundfilt[kmerhash] += 1
                else:
                    newabund = tuple([kmercount] + list(kmer.abund[1:]))
                    newkmer = KmerOfInterest(kmer.ksize, kmer.offset, newabund)
                    validated_kmers.append(newkmer)
                    self._valid[kmerhash] += 1
            record.annotations = validated_kmers
            if len(validated_kmers) == 0:
                self._novalidkmers_count += 1


def bam_filter_suppl(bam):
    for record in bam:
        if record.is_secondary or record.is_supplementary:
            continue
        yield record


def bam_paired_reader(bam):
    prev = None
    for record in bam_filter_suppl(bam):
        if not record.is_paired:
            assert prev is None
            yield record, None
            continue
        if prev is None:
            prev = record
            continue
        if prev.query_name == record.query_name:
            pr = prev.is_read1 and record.is_read2
            rp = prev.is_read2 and record.is_read1
            assert pr or rp
            if pr:
                yield prev, record
            else:
                yield record, prev
            prev = None
        else:
            yield prev, None
            if record.is_paired:
                prev = record
            else:
                yield record
                prev = None
