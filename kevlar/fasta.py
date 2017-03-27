#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import re
import screed


def parse(data):
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
    for defline, sequence in parse(data):
        seqid = defline[1:].replace('\t', ' ').split(' ')[0]
        assert seqid not in seqs
        seqs[seqid] = sequence
    return seqs


def parse_augmented_fastq(instream):
    record = None
    annot_kmers = dict()

    for line in instream:
        if line.startswith('@'):
            if record is not None:
                yield record, annot_kmers
                annot_kmers = dict()
            readid = line[1:].strip()
            seq = next(instream).strip()
            _ = next(instream)
            qual = next(instream).strip()
            record = screed.Record(name=readid, sequence=seq, quality=qual)
        elif line.endswith('#\n'):
            offset = len(line) - len(line.lstrip())
            line = line.strip()[:-1]
            abundances = re.split('\s+', line)
            kmer = abundances.pop(0)
            abundances = [int(a) for a in abundances]
            annot_kmers[offset] = (kmer, abundances)
    if record is not None:
        yield record, annot_kmers
