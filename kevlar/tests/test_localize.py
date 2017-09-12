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
from kevlar.localize import KevlarRefrSeqNotFound
from kevlar.localize import KevlarVariantLocalizationError
from kevlar.localize import KevlarNoReferenceMatchesError
from kevlar.localize import extract_regions
from kevlar.tests import data_file


def test_interval_set_simple():
    intervals = IntervalSet()
    assert intervals.get('chr1') == (None, None)

    intervals.add('chr1', 100, 25)
    intervals.add('chr1', 115, 25)
    assert intervals.get('chr1') == (100, 140)


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


def test_extract_region_basic():
    intervals = IntervalSet()
    intervals.add('bogus-genome-chr2', 10, 10)
    instream = open(data_file('bogus-genome/refr.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=0)]
    assert len(regions) == 1
    assert regions[0] == ('bogus-genome-chr2_10-20', 'GTTACATTAC')


def test_extract_region_basic_2():
    intervals = IntervalSet()
    intervals.add('simple', 49, 21)
    intervals.add('simple', 52, 21)
    intervals.add('simple', 59, 21)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=5)]
    assert len(regions) == 1
    assert regions[0] == ('simple_44-85',
                          'AATACTATGCCGATTTATTCTTACACAATTAAATTGCTAGT')


def test_extract_region_large_span():
    intervals = IntervalSet()
    intervals.add('simple', 100, 21)
    intervals.add('simple', 200, 21)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    with pytest.raises(KevlarVariantLocalizationError) as vle:
        _ = [r for r in extract_regions(instream, intervals, maxspan=100)]
    assert 'variant spans 221 bp (max 100)' in str(vle)

    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals)]
    assert len(regions) == 1
    assert regions[0][0] == 'simple_50-271'


def test_extract_region_missing_seq():
    intervals = IntervalSet()
    intervals.add('simple', 100, 21)
    intervals.add('simple', 200, 21)
    intervals.add('TheCakeIsALie', 42, 21)
    intervals.add('TheCakeIsALie', 100, 21)
    intervals.add('TheCakeIsALie', 77, 21)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    with pytest.raises(KevlarRefrSeqNotFound) as rnf:
        _ = [r for r in extract_regions(instream, intervals)]
    assert 'TheCakeIsALie' in str(rnf)


def test_extract_region_boundaries():
    intervals = IntervalSet()
    intervals.add('simple', 15, 31)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=20)]
    assert len(regions) == 1
    assert regions[0][0] == 'simple_0-66'

    intervals = IntervalSet()
    intervals.add('simple', 925, 31)
    intervals.add('simple', 955, 31)
    intervals.add('simple', 978, 31)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=20)]
    assert len(regions) == 1
    assert regions[0][0] == 'simple_905-1000'


def test_main(capsys):
    contig = data_file('localize-contig.fa')
    refr = data_file('localize-refr.fa')
    arglist = ['localize', '--ksize', '23', contig, refr]
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
