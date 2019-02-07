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
import kevlar
from kevlar.sequence import Record, KmerOfInterest
from kevlar.sequence import write_record, parse_augmented_fastx
from khmer import Counttable
from networkx import Graph, connected_components
from re import search


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
    partmatch = search(r'kvcc=(\d+)', readname)
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
