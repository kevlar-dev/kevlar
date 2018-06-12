# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

# distutils: language = c
# distutils: sources = src/align.c
# cython: c_string_type=str, c_string_encoding=ascii

import kevlar

cdef extern from 'align.h':
    void align(const char *target, const char *query, int match, int mismatch,
               int gapopen, int gapextend, char *cigar, int *score)


def contig_align(target, query, int match=1, int mismatch=2, int gapopen=5,
                 int gapextend=0):
    cdef char cigar[4096]
    cdef int score
    align(target, query, match, mismatch, gapopen, gapextend, cigar, &score)
    return cigar, score


def align_both_strands(target, query, int match=1, int mismatch=2,
                       int gapopen=5, int gapextend=0):
    cigar1, score1 = contig_align(
        target.sequence, query.sequence, match, mismatch, gapopen, gapextend
    )
    cigar2, score2 = contig_align(
        target.sequence, kevlar.revcom(query.sequence), match, mismatch,
        gapopen, gapextend
    )
    if score2 > score1:
        cigar = cigar2
        score = score2
        strand = -1
    else:
        cigar = cigar1
        score = score1
        strand = 1
    return score, cigar, strand
