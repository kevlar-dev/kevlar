#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import namedtuple
import sys
import kevlar


def ReadPair(object):
    def __init__(self, read1, read2, sharedkmer):
        self._r1 = read1
        self._r2 = read2
        self._seedkmer = sharedkmer
        self._k1 = None
        self._k2 = None
        self._merged = None

        self.head = None
        self.tail = None
        self.headkmer = None
        self.tailkmer = None
        self.overlap = None
        self.sameorient = None

        self.validate()

    @property
    def headseq(self):
        if self.sameorient:
            return self.head.sequence
        else:
            return kevlar.revcom(self.head.sequence)

    def __str__(self):
        return '{tailseq}\n{offset}{headseq}'.format(
            tailseq=self.tail.sequence, offset=' ' * self.offset,
            headseq=self.headseq
        )

    def incompatible(self):
        return self._k1 is None or self.overlap is None

    def check_kmer_freq(self):
        """Ignore any k-mers that occur multiple times in either read."""
        rckmer = kevlar.revcom(self._seedkmer)
        matches1 = [
            k for k in self._r1.annotations
            if kevlar.same_seq(self._r1.ikmerseq(k), self._seedkmer, rckmer)
        ]
        matches2 = [
            k for k in self._r2.annotations
            if kevlar.same_seq(self._r2.ikmerseq(k), self._seedkmer, rckmer)
        ]
        nmatches1 = len(matches1)
        nmatches2 = len(matches2)
        assert nmatches1 > 0 and nmatches1 > 0, (nmatches1, nmatches2)
        if nmatches1 == 1 and nmatches2 == 1:
            self._k1 = matches1[0]
            self._k2 = matches2[0]
            k1seq = self._r1.ikmerseq(self._k1)
            k2seq = self._r2.ikmerseq(self._k2)
            assert kevlar.same_seq(k1seq, k2seq)

    def set_head_and_tail(self):
        """Determine which read is `head` and which is `tail`.

        The `tail` read is situated to the left and always retains its
        orientation. The `head` read is situated to the right and will be
        reverse complemented if needed to match the `tail` read's orientation.

        The criteria for determining heads and tails is as follows.
        - the longer read is the tail
        - if the reads are of equal length, the read with the higher ikmer
          offset is the tail
        - if the reads are of equal length and have the same ikmer offset, the
          read whose name/ID is lexicographically larger is the tail
        """
        if len(self._r1) - len(self._r2) != 0:  # read length
            if len(self._r1) > len(self._r2):
                self.tail = self._r1
                self.head = self._r2
            else:
                self.tail = self._r2
                self.head = self._r1
        elif self._k1.offset != self._k2.offset:  # k-mer position
            if self._k1.offset > self._k2.offset:
                self.tail = self._r1
                self.head = self._r2
            else:
                self.tail = self._r2
                self.head = self._r1
        else:  # read names
            if self._r1.name > self._r2.name:
                self.tail = self._r1
                self.head = self._r2
            else:
                self.tail = self._r2
                self.head = self._r1
        if self.tail == self._r1:
            self.tailkmer = self._k1
            self.headkmer = self._k2
        else:
            self.tailkmer = self._k2
            self.headkmer = self._k1

    def calc_orientation_and_offset(self):
        tailkseq = self.tail.ikmerseq(self.tailkmer)
        headkseq = self.head.ikmerseq(self.headkmer)
        self.sameorient = tailkseq == headkseq
        tailpos = self.tailkmer.offset
        headpos = self.headkmer.offset
        ksize = self.headkmer.ksize
        if not self.sameorient:
            headpos = len(self.head) - (self.headkmer.offset + ksize)
        self.offset = headpos - tailpos

    @property
    def mergedseq(self):
        return self._merged

    def _merge(self):
        if self.headseq in self.tail.sequence:
            self._merged = pair.tail.sequence
        headindex = len(tailseq) - offset
        headsuffix = headseq[headindex:]
        tailprefix = tailseq[offset:offset+pair.overlap]
        assert tailprefix == headseq[:headindex], \
            'error: attempted to assemble incompatible reads'
        self._merged = pair.tail.sequennce + headsuffix

    def validate(self):
        self.check_kmer_freq()
        if self.incompatible:
            return
        self.set_head_and_tail()
        self.calc_orientation_and_offset()
        self._merge()
