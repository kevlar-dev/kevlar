# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

# distutils: language = c
# distutils: sources = src/align.c
# cython: c_string_type=str, c_string_encoding=ascii

cdef extern from 'align.h':
    void align(const char *target, const char *query, int match, int mismatch,
               int gapopen, int gapextend, char *cigar)

def contig_align(str target, str query, int match=1, int mismatch=2,
                 int gapopen=5, int gapextend=0):
    cdef char cigar[4096];
    align(target, query, match, mismatch, gapopen, gapextend, cigar)
    return cigar
