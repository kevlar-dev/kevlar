#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import os.path
import re
from subprocess import Popen, PIPE, check_call
import sys
from tempfile import TemporaryFile

import khmer
import kevlar
import pysam
import screed


class KevlarBWAError(RuntimeError):
    """Raised if a delegated BWA call fails for any reason."""
    pass


class KevlarInvalidCutoutDeflineError(ValueError):
    pass


class KevlarDeflineSequenceLengthMismatchError(RuntimeError):
    pass


def autoindex(refrfile, logstream=sys.stderr):
    if not os.path.isfile(refrfile):
        message = 'reference file {:s} does not exist'.format(refrfile)
        raise KevlarBWAError(message)

    bwtfile = refrfile + '.bwt'
    if os.path.isfile(bwtfile):
        return

    message = 'WARNING: BWA index not found for "{:s}"'.format(refrfile)
    message += ', indexing now'
    print('[kevlar::reference]', message, file=logstream)

    try:
        check_call(['bwa', 'index', refrfile])
    except Exception as err:  # pragma: no cover
        raise KevlarBWAError('Could not run "bwa index"') from err


def bwa_align(cmdargs, seqstring=None, seqfilename=None):
    if (not seqstring) is (not seqfilename):
        raise Exception('supply sequence string or file, not both')
    with TemporaryFile() as samfile:
        kmerseqs = dict()
        if seqstring:
            bwaproc = Popen(cmdargs, stdin=PIPE, stdout=samfile, stderr=PIPE,
                            universal_newlines=True)
            stdout, stderr = bwaproc.communicate(input=seqstring)
        else:
            bwaproc = Popen(cmdargs, stdout=samfile, stderr=PIPE)
            stdout, stderr = bwaproc.communicate()
        if bwaproc.returncode != 0:
            print(stderr, file=sys.stderr)
            raise KevlarBWAError('problem running BWA')
        samfile.seek(0)
        sam = pysam.AlignmentFile(samfile, 'r')
        for record in sam:
            if record.is_unmapped:
                continue
            seqid = sam.get_reference_name(record.reference_id)
            seq = record.seq
            if seq:
                kmerseqs[record.query_name] = seq
            else:
                seq = kmerseqs[record.query_name]
            yield seqid, record.reference_start, record.reference_end, seq


class ReferenceCutout(object):
    """An interval of the reference genome matched by a variant contig.

    Kevlar identifies, filters, and partitions reads containing novel k-mers,
    and then assembles these reads into contigs representing putative variants.
    This class is for handling small regions of the reference genome that have
    similarity to one of these contigs, to facilitate variant calling.

    Reference cutouts are operationally defined as follows.

    - decompose the contig into its constituent k-mers (seeds); note, the value
      of k used here is not necessarily the same k that is used for variant
      discovery
    - find all perfect reference genome matches of all seeds
    - sort each seed match by genomic position
    - split two adjacent matches into distinct bins if they are on different
      reference sequences (chromosomes) or if they are separated by more than X
      nucleotides
    - compute the smallest interval spanning all matches in a bin
    - extend the interval by delta nucleotides and designate as a reference
      cutout
    """
    def __init__(self, defline=None, sequence=None):
        self.defline = defline
        self.sequence = sequence
        self._seqid = None
        self._startpos = None
        self._endpos = None
        if defline:
            self.parse_defline(defline)

    def __len__(self):
        return self._endpos - self._startpos

    def parse_defline(self, defline):
        match = re.search(r'(\S+)_(\d+)-(\d+)', defline)
        if not match:
            raise KevlarInvalidCutoutDeflineError(defline)
        self._seqid = match.group(1)
        self._startpos = int(match.group(2))
        self._endpos = int(match.group(3))
        if not self.sequence:
            return
        if len(self) != len(self.sequence):
            message = 'defline length: {:d}, sequence length: {:d}'.format(
                len(self), len(self.sequence)
            )
            raise KevlarDeflineSequenceLengthMismatchError(message)

    @property
    def interval(self):
        return self._seqid, self._startpos, self._endpos

    def local_to_global(self, coordinate):
        return self._startpos + coordinate


def load_refr_cutouts(instream):
    for defline, sequence in kevlar.seqio.parse_fasta(instream):
        yield ReferenceCutout(defline[1:], sequence)
