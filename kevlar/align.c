#include <Python.h>
#include "ksw2.h"

// Stolen shamelessly from ksw2/cli.c
unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// Stolen shamelessly from ksw2/cli.c
static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}

static PyObject *align_call(PyObject *self, PyObject *args)
{
    /**
     * Align `query` to `target` using the `ksw_extz` algorithm. See
     * https://github.com/lh3/ksw2 for more info.
     *
     * @param target    reference genome sequence
     * @param query     contig assembled from variant-associated reads
     * @returns         CIGAR string of the alignment
     */
    const char *target;
    const char *query;
    char cigar[512];
    size_t ci = 0;

	int8_t score_match = 1;
    int8_t penalty_mismatch = 2;
    int8_t penalty_gapopen = 5
    int8_t penalty_gap_extend = 0.5;
    int8_t alphabet_size = 5;

    int8_t matrix[25];
    ksw_extz_t ez;

    if (!PyArg_ParseTuple(args, "ss", &target, &query)) {
        return NULL;
    }

    memset(&ez, 0, sizeof(ksw_extz_t));
    ksw_gen_simple_mat(5, matrix, score_match, -penalty_mismatch);

    ksw_extz(
        NULL, // memory pool
        strlen(query), query, strlen(target), target, alphabet_size, matrix,
        penalty_gapopen, penalty_gap_extend,
        -1, // bandwidth
        -1, // zdrop
        0, //flag
        ez
    );

    // Stolen shamelessly from ksw2/cli.c
    for (int i = 0; i < ez->n_cigar; ++i) {
		ci += sprintf(
            cigar + ci, "%d%c", ez->cigar[i] >> 4,"MID"[ez->cigar[i]&0xf]
        );
    }

    return PyUnicodeFromString(cigar);
}
