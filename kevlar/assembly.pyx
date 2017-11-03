# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

# distutils: language = c
# distutils: sources = src/assemble.c
# cython: c_string_type=str, c_string_encoding=ascii
from libc.stdint cimport int32_t
from libc.stdlib cimport malloc, free
from libc.string cimport strdup

cdef extern from 'fml.h':
    ctypedef struct bseq1_t:
        int32_t l_seq
        char *seq
        char *qual

    ctypedef struct fml_opt_t:
        pass

    ctypedef struct fml_utg_t:
        char *seq

    cdef void fml_opt_init(fml_opt_t *)
    cdef fml_utg_t *fml_assemble(const fml_opt_t *, int, bseq1_t *, int *)
    cdef void fml_utg_destroy(int, fml_utg_t *)

def fml_asm(records):
    cdef int num_unitigs
    cdef int numreads = len(records)
    cdef bseq1_t *reads = <bseq1_t *>malloc(sizeof(bseq1_t) * numreads)
    for i, read in enumerate(records):
        reads[i].seq = strdup(read.sequence)
        reads[i].l_seq = len(read.sequence)
        reads[i].qual = NULL

    cdef fml_opt_t opt
    fml_opt_init(&opt)
    cdef fml_utg_t *unitigs = fml_assemble(&opt, numreads, reads, &num_unitigs)
    for i in range(num_unitigs):
        yield unitigs[i].seq
    fml_utg_destroy(num_unitigs, unitigs)
