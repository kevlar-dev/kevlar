# -----------------------------------------------------------------------------
# Copyright (c) 2019 Battelle National Biodefense Institute
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

# distutils: language = c++
# distutils: sources = src/signature.cpp
# cython: c_string_type=str, c_string_encoding=ascii

from libc.stdint cimport uint8_t
from libcpp.vector cimport vector
import khmer
import kevlar

cdef extern from 'signature.h':
    cdef cppclass NovelKmer:
        unsigned ksize
        unsigned offset
        vector[uint8_t] abunds

    cdef cppclass NovelRead:
        oxli::read_parsers:Read read
        vector[NovelKmer] annotations

    int next_novel_read(NovelRead buffer, Parser &parser,
                        oxli::hashtable::Counttable &casecounts,
                        std::vector<oxli::hashtable::Counttable> &controlcounts,
                        unsigned casemin, unsigned ctrlmax)


def get_novel_reads(parser, casecounts, controlcounts, casemin=5, ctrlmax=1):
    cdef NovelRead r
    while True:
        result = next_novel_read(r, parser, casecounts, controlcounts, casemin,
                                 ctrlmax)
        if result == 0:
            return
        yield kevlar.sequence.Record(
            r.read.name, r.read.sequence, r.read.quality,
            annotations=r.annotations
        )
