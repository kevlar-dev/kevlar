#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import namedtuple
import sys
import kevlar


class ReadWithKmer(object):
    def __init__(self, read, kmerseq):
        self.read = read
        self.kmerocc = list()
        kmerseqrc = kevlar.revcom(kmerseq)
        for ikmer in read.annotations:
            if kevlar.same_seq(self.read.ikmerseq(ikmer), kmerseq, kmerseqrc):
                self.kmerocc.append(ikmer)
        self.kmer = self.kmerocc[0] if len(self.kmerocc) == 1 else None

    def __len__(self):
        return len(self.read)

    @property
    def num_occurrences(self):
        return len(self.kmerocc)

    @property
    def kmerseq(self):
        if self.kmer is None:
            return None
        return self.read.ikmerseq(self.kmer)

    @property
    def offset(self):
        return self.kmer.offset

    @property
    def name(self):
        return self.read.name

    @property
    def sequence(self):
        return self.read.sequence


class ReadPair(object):
    """Class for handling overlapping pairs of reads.

    Each read pair is anchored by 1 share k-mer at a time. These shared k-mers
    are used to compute relative orientation, offset, sequence compatibility/
    identity, and to "assemble" (merge) the two sequences.
    """
    def __init__(self, read1, read2, sharedkmer):
        self._r1 = ReadWithKmer(read1, sharedkmer)
        self._r2 = ReadWithKmer(read2, sharedkmer)
        self._r1rc = ReadWithKmer(read1.revcom(), sharedkmer)
        self._r2rc = ReadWithKmer(read2.revcom(), sharedkmer)
        self._seedkmer = sharedkmer
        self._merged = None

        self.head = None
        self.tail = None
        self.headkmer = None
        self.tailkmer = None
        self.overlap = None
        self.sameorient = None

        self.validate()

    def __str__(self):
        return '{tailseq}\n{koffset}{match}\n{offset}{headseq}'.format(
            tailseq=self.tail.read.sequence, koffset=' ' * self.tail.offset,
            match='|' * self.tail.kmer.ksize, offset=' ' * self.offset,
            headseq=self.head.read.sequence
        )

    @property
    def incompatible(self):
        return self._merged is None

    def assign_by_largest_kmer_offset(self):
        """Assign head and tail based on largest k-mer offset."""
        # Enumerate all possible arragements of the read pair based on the
        # relative orientation of the reads.
        if self.sameorient:
            arrangements = [(self._r1, self._r2), (self._r1rc, self._r2rc)]
        else:
            arrangements = [(self._r1, self._r2rc), (self._r1rc, self._r2)]

        # If both arrangements have the same largest k-mer offset, no assigment
        # can be made.
        offsets = [read.kmer.offset for arr in arrangements for read in arr]
        if len(set(offsets)) == 1:
            return

        # Determine the optimal arrangement and assign head and tail.
        if offsets[0] > offsets[1]:
            optimal_arrangement = arrangements[0]
        else:
            optimal_arrangement = arrangements[1]
        self.tail = max(optimal_arrangement, key=lambda r: r.kmer.offset)
        self.head = min(optimal_arrangement, key=lambda r: r.kmer.offset)

    def assign_by_read_length(self):
        if len(self._r1) == len(self._r2):
            return
        elif len(self._r1) > len(self._r2):
            self.tail = self._r1
            self.head = self._r2 if self.sameorient else self._r2rc
        else:
            self.tail = self._r2 if self.sameorient else self._r2rc
            self.head = self._r1

    def assign_by_read_name(self):
        if len(self._r1.read.name) < len(self._r2.read.name):
            self.tail = self._r1
            self.head = self._r2 if self.sameorient else self._r2rc
        else:
            self.tail = self._r2 if self.sameorient else self._r2rc
            self.head = self._r1

    def set_head_and_tail(self):
        """Determine which read is `head` and which is `tail`.

        The `tail` read is situated to the left and always retains its
        orientation. The `head` read is situated to the right and will be
        reverse complemented if needed to match the `tail` read's orientation.

        The criteria for determining heads and tails is as follows.
        - the read with the higher k-mer offset is the tail; all possible
          arrangements are considered
        - if all arrangements have the same k-mer offset, the longer read is
          the tail
        - if the reads are of equal length, the read whose name/ID is
          lexicographically smaller is the tail
        """
        self.assign_by_largest_kmer_offset()
        if self.tail is None:
            self.assign_by_read_length()
        if self.tail is None:
            self.assign_by_read_name()
        assert self.tail is not None

    def calc_offset(self):
        self.offset = self.tail.offset - self.head.offset
        self.overlap = len(self.tail) - self.offset

    @property
    def mergedseq(self):
        return self._merged

    def _merge(self):
        tailseq = self.tail.read.sequence
        headseq = self.head.read.sequence
        if headseq in tailseq:
            self._merged = tailseq
            return
        elif tailseq in headseq:
            self._merged = tailseq
            return

        headindex = len(tailseq) - self.offset
        headsuffix = headseq[headindex:]
        tailprefix = tailseq[self.offset:self.offset+self.overlap]
        if tailprefix == headseq[:headindex]:
            self._merged = tailseq + headsuffix

    def validate(self):
        if self._r1.num_occurrences != 1 or self._r2.num_occurrences != 1:
            return
        self.sameorient = self._r1.kmerseq == self._r2.kmerseq
        self.set_head_and_tail()
        self.calc_offset()
        self._merge()
