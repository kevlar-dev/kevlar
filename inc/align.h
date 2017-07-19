//------------------------------------------------------------------------------
//Copyright (c) 2017 The Regents of the University of California
//
// This file is part of kevlar (http://github.com/standage/kevlar) and is
// licensed under the MIT license: see LICENSE.
//------------------------------------------------------------------------------

/**
 * Align `query` to `target` using the `ksw_extz` algorithm. See
 * https://github.com/lh3/ksw2 for more info.
 *
 * @param target      reference genome sequence
 * @param query       contig assembled from variant-associated reads
 * @param match       score for a nucleotide match
 * @param mismatch    penalty for a nucleotide mismatch
 * @param gapopen     gap open penalty
 * @param gapextend   gap extension penalty
 * @returns           CIGAR string of the alignment; calling function is
                      responsible for freeing heap-allocated memory
 */
char *align(const char *target, const char *query, int match,
            int mismatch, int gapopen, int gapextend);
