#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from io import StringIO
import pytest
import screed
import kevlar
from kevlar.tests import data_file


@pytest.mark.parametrize('infile', [
    ('smallseq.fa.gz'),
    ('smallseq.augfasta'),
    ('smallseq.fq'),
])
def test_get_unique_kmers(infile):
    infile = data_file(infile)
    kmers = set([k for k in kevlar.localize.get_unique_kmers(infile, ksize=9)])
    testkmers = set(
        ['TTAATTGGC', 'CTTAATTGG', 'TAATTGGCC', 'ATTACCGGT',
         'TTACCGGTA', 'CCTTAATTG', 'GCCTTAATT', 'GGCCTTAAT']
    )
    print(sorted(kmers))
    print(sorted(testkmers))
    assert kmers == testkmers


def test_unique_kmer_string():
    infile = data_file('smallseq.augfasta')
    fastastring = kevlar.localize.unique_kmer_string(infile, ksize=9)
    fastafile = StringIO(fastastring)
    kmers = set([s for d, s in kevlar.seqio.parse_fasta(fastafile)])
    testkmers = set(
        ['TTAATTGGC', 'CTTAATTGG', 'TAATTGGCC', 'ATTACCGGT',
         'TTACCGGTA', 'CCTTAATTG', 'GCCTTAATT', 'GGCCTTAAT']
    )
    print(sorted(kmers))
    print(sorted(testkmers))
    assert kmers == testkmers


def test_select_region():
    # Different sequences
    matches = [('chr1', 100), ('chr2', 450)]
    assert kevlar.localize.select_region(matches) is None

    # Too distant
    matches = [('chr1', 100), ('chr1', 45000)]
    assert kevlar.localize.select_region(matches, maxdiff=1000) is None

    # On the same sequence, close together, passes!
    matches = [('chr1', 4000), ('chr1', 4032), ('chr1', 3990)]
    region = ('chr1', 3890, 4133)
    assert kevlar.localize.select_region(matches, delta=100) == region
