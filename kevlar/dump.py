#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import pysam
import kevlar
from kevlar.seqio import bam_paired_reader
from kevlar.sequence import Record, write_record


def perfectmatch(record, bam, refrseqs):
    """
    Determine if the alignment record is a perfect match against the reference.

    The `record` is a pysam alignment object, `bam` is a pysam BAM parser
    object, and `refrseqs` is a dictionary of sequences indexed by their
    sequence IDs.
    """
    matchcigar = '{:d}M'.format(record.rlen)
    if record.cigarstring == matchcigar:
        seq = refrseqs[bam.get_reference_name(record.tid)]
        refrsubseq = str(seq[record.pos:record.pos+record.rlen])
        if refrsubseq.upper() == record.seq.upper():
            return True
    return False


def readname(record):
    """Create a Fastq read name, using suffixes for paired reads as needed."""
    name = record.query_name
    if record.flag & 1:
        # Logical XOR: if the read is paired, it should be first in pair
        # or second in pair, not both.
        assert (record.flag & 64) != (record.flag & 128)
        suffix = '/1' if record.flag & 64 else '/2'
        if not name.endswith(suffix):
            name += suffix
    return name


def keepers(record1, record2, bam, refrseqs=None, pairmode='split'):
    assert pairmode in ('split', 'keep', 'drop'), pairmode
    if record2 is None:  # single end behavior
        keep = not refrseqs or not perfectmatch(record1, bam, refrseqs)
        if keep:
            return [record1]
        else:
            return []

    r1keep = not refrseqs or not perfectmatch(record1, bam, refrseqs)
    r2keep = not refrseqs or not perfectmatch(record2, bam, refrseqs)
    if r1keep == r2keep:
        if r1keep and r2keep:
            return [record1, record2]
        else:
            return []
    else:
        if pairmode == 'split':
            if r1keep:
                return [record1]
            else:
                assert r2keep
                return [record2]
        elif pairmode == 'keep':
            return [record1, record2]
        else:
            assert pairmode == 'drop', pairmode
            return []


def dump(bamstream, refrseqs=None, pairmode='split', upint=50000,
         logstream=sys.stderr):
    """
    Parse read alignments in BAM/SAM format.

    - bamstream: open file handle to the BAM/SAM file input
    - refrseqs: dictionary of reference sequences, indexed by sequence ID; if
      provided, perfect matches to the reference sequence will be discarded
    - strict: only keep paired end if it also lacks a perfect match to the
      reference genome
    - upint: update interval for progress indicator
    - logstream: file handle do which progress indicator will write output
    """
    bam = pysam.AlignmentFile(bamstream, 'rb')
    reader = bam_paired_reader(bam)
    for i, (record1, record2) in enumerate(reader, 1):
        if i % upint == 0:  # pragma: no cover
            print('...processed', i, 'pairs of records', file=logstream)
        for record in keepers(record1, record2, bam, refrseqs, pairmode):
            yield Record(
                name=readname(record), sequence=record.seq, quality=record.qual
            )


def main(args):
    fastq = kevlar.open(args.out, 'w')
    refr = None
    if args.refr:
        print('[kevlar::dump] Loading reference sequence', file=args.logfile)
        refrstream = kevlar.open(args.refr, 'r')
        refr = kevlar.seqio.parse_seq_dict(refrstream)
    for read in dump(args.reads, refr, args.pair_mode, logstream=args.logfile):
        write_record(read, fastq)
