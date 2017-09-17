#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import pytest
import screed
import kevlar
from kevlar.localize import IntervalSet
from kevlar.localize import (KevlarRefrSeqNotFoundError,
                             KevlarVariantLocalizationError,
                             KevlarNoReferenceMatchesError)
from kevlar.localize import (extract_regions, get_unique_kmers,
                             unique_kmer_string)
from kevlar.tests import data_file


def test_interval_set_simple():
    intervals = IntervalSet()
    assert intervals.get_spans('chr1') is None

    intervals.add('chr1', 100, 125)
    intervals.add('chr1', 115, 140)
    assert intervals.get_spans('chr1') == [(100, 140)]

    intervals.add('chr2', 200, 225)
    intervals.add('chr2', 205, 230)
    intervals.add('chr2', 207, 232)
    intervals.add('chr2', 235008, 235033)
    intervals.add('chr2', 235075, 235100)
    assert intervals.get_spans('chr2') == [(200, 232), (235008, 235100)]


@pytest.mark.parametrize('infile', [
    ('smallseq.fa.gz'),
    ('smallseq.augfasta'),
    ('smallseq.fq'),
])
def test_get_unique_kmers(infile):
    infile = data_file(infile)
    instream = kevlar.parse_augmented_fastx(kevlar.open(infile, 'r'))
    kmers = set([k for k in get_unique_kmers(instream, ksize=9)])
    testkmers = set(
        ['TTAATTGGC', 'CTTAATTGG', 'TAATTGGCC', 'ATTACCGGT',
         'TTACCGGTA', 'CCTTAATTG', 'GCCTTAATT', 'GGCCTTAAT']
    )
    print(sorted(kmers))
    print(sorted(testkmers))
    assert kmers == testkmers


def test_unique_kmer_string():
    infile = data_file('smallseq.augfasta')
    instream = kevlar.parse_augmented_fastx(kevlar.open(infile, 'r'))
    fastastring = unique_kmer_string(instream, ksize=9)
    fastafile = StringIO(fastastring)
    kmers = set([s for d, s in kevlar.seqio.parse_fasta(fastafile)])
    testkmers = set(
        ['TTAATTGGC', 'CTTAATTGG', 'TAATTGGCC', 'ATTACCGGT',
         'TTACCGGTA', 'CCTTAATTG', 'GCCTTAATT', 'GGCCTTAAT']
    )
    print(sorted(kmers))
    print(sorted(testkmers))
    assert kmers == testkmers


def test_extract_regions_basic():
    intervals = IntervalSet()
    intervals.add('bogus-genome-chr2', 10, 20)
    instream = open(data_file('bogus-genome/refr.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=0)]
    assert len(regions) == 1
    assert regions[0] == ('bogus-genome-chr2_10-20', 'GTTACATTAC')


def test_extract_regions_basic_2():
    intervals = IntervalSet()
    intervals.add('simple', 49, 49 + 21)
    intervals.add('simple', 52, 52 + 21)
    intervals.add('simple', 59, 59 + 21)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=5)]
    assert len(regions) == 1
    assert regions[0] == ('simple_44-85',
                          'AATACTATGCCGATTTATTCTTACACAATTAAATTGCTAGT')


def test_extract_regions_large_span():
    intervals = IntervalSet()
    intervals.add('simple', 100, 121)
    intervals.add('simple', 200, 221)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    seqids = [r[0] for r in extract_regions(instream, intervals, maxdiff=50)]
    assert seqids == ['simple_75-146', 'simple_175-246']

    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=50)]
    assert len(regions) == 1
    assert regions[0][0] == 'simple_50-271'


def test_extract_regions_missing_seq():
    intervals = IntervalSet()
    intervals.add('simple', 100, 121)
    intervals.add('simple', 200, 221)
    intervals.add('TheCakeIsALie', 42, 63)
    intervals.add('TheCakeIsALie', 100, 121)
    intervals.add('TheCakeIsALie', 77, 98)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    with pytest.raises(KevlarRefrSeqNotFoundError) as rnf:
        _ = [r for r in extract_regions(instream, intervals)]
    assert 'TheCakeIsALie' in str(rnf)


def test_extract_regions_boundaries():
    intervals = IntervalSet()
    intervals.add('simple', 15, 46)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=20)]
    assert len(regions) == 1
    assert regions[0][0] == 'simple_0-66'

    intervals = IntervalSet()
    intervals.add('simple', 925, 925 + 31)
    intervals.add('simple', 955, 955 + 31)
    intervals.add('simple', 978, 978 + 31)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=20)]
    assert len(regions) == 1
    assert regions[0][0] == 'simple_905-1000'


def test_main(capsys):
    contig = data_file('localize-contig.fa')
    refr = data_file('localize-refr.fa')
    arglist = ['localize', '--ksize', '23', '--delta', '50', contig, refr]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.localize.main(args)
    out, err = capsys.readouterr()
    assert '>seq1_10-191' in out


def test_main_no_matches():
    contig = data_file('localize-contig-bad.fa')
    refr = data_file('localize-refr.fa')
    arglist = ['localize', '--ksize', '23', contig, refr]
    args = kevlar.cli.parser().parse_args(arglist)
    with pytest.raises(KevlarNoReferenceMatchesError) as nrm:
        kevlar.localize.main(args)
