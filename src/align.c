//------------------------------------------------------------------------------
//Copyright (c) 2017 The Regents of the University of California
//
// This file is part of kevlar (http://github.com/standage/kevlar) and is
// licensed under the MIT license: see LICENSE.
//------------------------------------------------------------------------------

#include <string.h>
#include "ksw2.h"

/**
 * Encode the given nucleotide sequence.
 *
 * @param sequence    nucleotide sequence
 * @returns           pointer to an array of uint8_t; array is allocated on the
                      heap, and calling function is responsible to free the
                      memory
 */
static uint8_t* encode(const char *sequence)
{
    uint8_t table[256];
    memset(table, 4, 256);
    table[(int)'A'] = table[(int)'a'] = 0;
    table[(int)'C'] = table[(int)'c'] = 1;
    table[(int)'G'] = table[(int)'g'] = 2;
    table[(int)'T'] = table[(int)'t'] = 3;

    uint8_t *int_seq = (uint8_t *)malloc(strlen(sequence));
    for (size_t i = 0; i < strlen(sequence); i++) {
        int_seq[i] = table[(uint8_t)sequence[i]];
    }

    return int_seq;
}


// See align.h for documentation.
char* align(const char *target, const char *query, int match,
            int mismatch, int gapopen, int gapextend)
{
    int8_t a = match; // a > 0
    int8_t b = mismatch < 0 ? mismatch : -mismatch; // b < 0
    int8_t matrix[25] = {
        a,b,b,b,0,
        b,a,b,b,0,
        b,b,a,b,0,
        b,b,b,a,0,
        0,0,0,0,0
    };

    uint8_t *target_enc = encode(target);
    uint8_t *query_enc = encode(query);

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    ksw_extz(
        NULL, // memory pool
        strlen(query), query_enc,
        strlen(target), target_enc,
        5, // alphabet size, |{A,C,G,T,N}|
        matrix, gapopen, gapextend,
        -1, // bandwidth
        -1, // zdrop
        0, //flag
        &ez
    );

    // Stolen shamelessly from ksw2/cli.c
    char cigar[4096];
    size_t ci = 0;
    for (size_t i = 0; i < ez.n_cigar; ++i) {
        ci += sprintf(
            cigar + ci, "%d%c", ez.cigar[i] >> 4,"MID"[ez.cigar[i]&0xf]
        );
    }

    char *cigarstring = malloc(strlen(cigar) + 1);
    strcpy(cigarstring, cigar);
    free(query_enc);
    free(target_enc);
    return cigarstring;
}
