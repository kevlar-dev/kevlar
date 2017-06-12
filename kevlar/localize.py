#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from subprocess import Popen, PIPE
from tempfile import TemporaryFile
import sys

import kevlar
import khmer
import pysam


class KevlarBWAError(RuntimeError):
    pass


class KevlarNoReferenceMatchesError(ValueError):
    pass


class KevlarVariantLocalizationError(ValueError):
    pass


class KevlarRefrSeqNotFound(ValueError):
    pass


def get_unique_kmers(infile, ksize=31):
    ct = khmer._Counttable(ksize, [1])
    kmers = set()
    instream = open(infile, 'r')
    for defline, sequence in kevlar.seqio.parse_fasta(instream):
        for kmer in ct.get_kmers(sequence):
            minkmer = kevlar.revcommin(kmer)
            if minkmer not in kmers:
                kmers.add(minkmer)
                yield kmer


def unique_kmer_string(infile, ksize=31):
    output = ''
    for n, kmer in enumerate(get_unique_kmers(infile, ksize)):
        output += '>kmer{:d}\n{:s}\n'.format(n, kmer)
    return output


def get_exact_matches(infile, bwaindexfile, ksize=31):
    kmers = unique_kmer_string(infile, ksize)
    cmd = 'bwa mem -k {k} -T {k} {idx} -'.format(k=ksize, idx=bwaindexfile)
    cmdargs = cmd.split(' ')
    with TemporaryFile() as samfile:
        bwaproc = Popen(cmdargs, stdin=PIPE, stdout=samfile, stderr=PIPE,
                        universal_newlines=True)
        stdout, stderr = bwaproc.communicate(input=kmers)
        if bwaproc.returncode != 0:
            print(stderr, file=sys.stderr)
            raise KevlarBWAError('problem running BWA')
        samfile.seek(0)
        sam = pysam.AlignmentFile(samfile, 'r')
        for record in sam:
            if record.is_unmapped:
                continue
            seqid = sam.get_reference_name(record.reference_id)
            yield seqid, record.pos


def select_region(matchlist, maxdiff=1000, delta=100):
    seqids = set([s for s, p in matchlist])
    if len(seqids) > 1:
        return None

    minpos = min([p for s, p in matchlist])
    maxpos = max([p for s, p in matchlist])
    if maxpos - minpos > maxdiff:
        return None
    return seqids.pop(), minpos-delta, maxpos+delta+1


def extract_region(refr, seqid, start, end):
    for defline, sequence in kevlar.seqio.parse_fasta(refr):
        testseqid = defline[1:].split()[0]
        if seqid == testseqid:
            subseqid = '{}_{}-{}'.format(seqid, start, end)
            subseq = sequence[start:end]
            return subseqid, subseq
    raise KevlarRefrSeqNotFound()


def main(args):
    output = args.out if args.out else sys.stdout
    matchgen = get_exact_matches(args.contigs, args.refr, args.ksize)
    kmer_matches = [m for m in matchgen]
    if len(kmer_matches) == 0:
        raise KevlarNoReferenceMatchesError()

    region = select_region(kmer_matches, args.max_diff, args.delta)
    if region is None:
        raise KevlarVariantLocalizationError()
    else:
        seqid, start, end = region
        instream = kevlar.open(args.refr, 'r')
        subseqid, subseq = extract_region(instream, seqid, start, end)
        print('>', subseqid, '\n', subseq, sep='', file=output)
