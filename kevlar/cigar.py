#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import namedtuple
import kevlar
import re
import sys


AlignmentBlock = namedtuple('AlignmentBlock', 'length type target query')


class AlignmentTokenizer(object):
    def __init__(self, queryseq, targetseq, cigar):
        self._query = queryseq
        self._target = targetseq
        self._origcigar = cigar
        self._cigar = cigar

        self.blocks = list(self._tokenize())
        self._endcheck()

    def _tokenize(self):
        target = self._target
        query = self._query
        blocks = re.finditer(r'(\d+)([DIM])', self._origcigar)
        for block in blocks:
            length = int(block.group(1))
            blocktype = block.group(2)
            tseq, qseq = None, None
            if blocktype in ('M', 'D'):
                tseq = target[:length]
                target = target[length:]
            if blocktype in ('M', 'I'):
                qseq = query[:length]
                query = query[length:]
            yield AlignmentBlock(length, blocktype, tseq, qseq)
        assert target == ''
        assert query == ''

    def _endcheck(self):
        if len(self.blocks) < 3:
            return
        if self.blocks[-1].type != 'M' or self.blocks[-3].type != 'M':
            return

        if self.blocks[-2].type == 'D':
            prevseq = self.blocks[-2].target
            lastseq = self.blocks[-1].target
            endseq = self.blocks[-1].query
        else:
            prevseq = self.blocks[-2].query
            lastseq = self.blocks[-1].query
            endseq = self.blocks[-1].target
        longseq = prevseq + lastseq
        if longseq.startswith(endseq):
            self.blocks[-3] = AlignmentBlock(
                self.blocks[-3].length + self.blocks[-1].length, 'M',
                self.blocks[-3].target + self.blocks[-1].target,
                self.blocks[-3].query + self.blocks[-1].query,
            )
            del self.blocks[-1]
            newcigar = ''
            for block in self.blocks:
                newcigar += '{:d}{:s}'.format(block.length, block.type)
            self._cigar = newcigar
