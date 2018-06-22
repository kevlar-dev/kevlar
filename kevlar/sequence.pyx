# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

# cython: c_string_type=str, c_string_encoding=ascii

from collections import namedtuple
import kevlar

KmerOfInterest = namedtuple('KmerOfInterest', 'ksize offset abund')

def tostr(stringlike):
    try:
        stringlike = stringlike.decode('utf-8')
    except AttributeError:
        pass
    return stringlike


cdef class Record:
    cdef public str name
    cdef public str sequence
    cdef public str quality
    cdef public list annotations
    cdef public list mates

    def __init__(self, str name, str sequence, str quality=None,
                 list annotations=None, list mates=None):
        self.name = name
        self.sequence = sequence
        self.quality = quality
        self.annotations = list() if annotations is None else annotations
        self.mates = list() if mates is None else mates

    def __len__(self):
        return len(self.sequence)

    def add_mate(self, str mateseq):
        self.mates.append(mateseq)

    def annotate(self, str sequence, int offset, tuple abundances):
        checkseq = self.sequence[offset:offset+len(sequence)]
        assert checkseq == sequence, (checkseq, sequence)
        ikmer = KmerOfInterest(len(sequence), offset, abundances)
        self.annotations.append(ikmer)

    @property
    def id(self):
        return self.name.split()[0]

    def ikmerseq(self, ikmer):
        return self.sequence[ikmer.offset:ikmer.offset+ikmer.ksize]


def copy_record(record):
    qual = None
    if hasattr(record, 'quality') and record.quality is not None:
        qual = record.quality
    return Record(record.name, record.sequence, qual)


def write_record(Record record, outstream):
    if record.quality is not None:
        recstr = '@{name}\n{sequence}\n+\n{quality}\n'.format(
            name=record.name,
            sequence=record.sequence,
            quality=record.quality
        )
    else:
        recstr = '>{name}\n{sequence}\n'.format(
            name=record.name,
            sequence=record.sequence
        )
    if len(record.annotations) > 0:
        annstrs = list()
        for kmer in sorted(record.annotations, key=lambda k: k.offset):
            abundstr = ' '.join([str(a) for a in kmer.abund])
            annstr = '{padding}{seq}{margin}{abund}#'.format(
                padding=' ' * kmer.offset,
                seq=record.sequence[kmer.offset:kmer.offset+kmer.ksize],
                margin=' ' * 10,
                abund=abundstr,
            )
            annstrs.append(annstr)
        recstr += '\n'.join(annstrs) + '\n'
    if len(record.mates) > 0:
        matestrs = list()
        for mateseq in record.mates:
            matestr = '#mateseq={:s}#'.format(mateseq)
            matestrs.append(matestr)
        recstr += '\n'.join(matestrs) + '\n'
    try:
        outstream.write(bytes(recstr, 'ascii'))
    except TypeError:
        outstream.write(recstr)
