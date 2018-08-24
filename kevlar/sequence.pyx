# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

# cython: c_string_type=str, c_string_encoding=ascii

from collections import namedtuple
import re
import kevlar

KmerOfInterest = namedtuple('KmerOfInterest', 'ksize offset abund')

revcomptab = str.maketrans(
    'ATUGCYRSWKMBDHVNatugcyrswkmbdhvn',
    'TAACGRYSWMKVHDBNTAACGRYSWMKVHDBN'
)


def revcom(str sequence):
    return sequence.translate(revcomptab)[::-1]


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
    cdef public dict ikmers

    def __init__(self, str name, str sequence, str quality=None,
                 list annotations=None, list mates=None, dict ikmers=None):
        self.name = name
        self.sequence = sequence
        self.quality = quality
        self.mates = list() if mates is None else mates
        self.ikmers = dict()
        if annotations is None:
            self.annotations = list()
            self.ikmers = dict()
        else:
            self.annotations = annotations
            if ikmers is None:
                for kmer in annotations:
                    kmerseq = self.ikmerseq(kmer)
                    kmerseqrc = revcom(kmerseq)
                    self.ikmers[kmerseq] = kmer
                    self.ikmers[kmerseqrc] = kmer
            else:
                self.ikmers = ikmers


    def __len__(self):
        return len(self.sequence)

    def add_mate(self, str mateseq):
        self.mates.append(mateseq)

    def annotate(self, str sequence, int offset, tuple abundances):
        checkseq = self.sequence[offset:offset+len(sequence)]
        assert checkseq == sequence, (checkseq, sequence)
        ikmer = KmerOfInterest(len(sequence), offset, abundances)
        self.annotations.append(ikmer)
        self.ikmers[sequence] = ikmer
        self.ikmers[revcom(sequence)] = ikmer

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


def print_augmented_fastx(Record record, outstream):
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


cpdef write_record(Record record, outstream):
    print_augmented_fastx(record, outstream)


def parse_augmented_fastx(instream):
    """Read augmented Fast[q|a] records into memory.

    The parsed records will have .name, .sequence, and .quality defined (unless
    it's augmented Fasta), as well as a list of interesting k-mers. See
    http://kevlar.readthedocs.io/en/latest/formats.html#augmented-sequences for
    more information.
    """
    cdef Record record = None
    cdef str readname
    cdef str seq
    cdef str qual
    cdef str mateseq
    cdef str firstchar
    cdef int offset
    cdef str kmer

    for line in instream:
        if line.strip() == '':
            continue
        firstchar = line[0]
        if firstchar in ('@', '>'):
            if record is not None:
                yield record
            readname = line[1:].strip()
            seq = next(instream).strip()
            if firstchar == '@':
                _ = next(instream)
                qual = next(instream).strip()
            else:
                qual = None
            record = Record(name=readname, sequence=seq, quality=qual)
        elif line.endswith('#\n'):
            if line.startswith('#mateseq='):
                mateseq = re.search(r'^#mateseq=(\S+)#\n$', line).group(1)
                record.add_mate(mateseq)
                continue
            offset = len(line) - len(line.lstrip())
            line = line.strip()[:-1]
            abundlist = re.split(r'\s+', line)
            kmer = abundlist.pop(0)
            abundances = tuple([int(a) for a in abundlist])
            record.annotate(kmer, offset, abundances)
        else:
            raise Exception(line)
    yield record
