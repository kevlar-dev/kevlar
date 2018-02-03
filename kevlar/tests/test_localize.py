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
from kevlar.localize import SeedMatchSet
from kevlar.localize import KevlarRefrSeqNotFoundError
from kevlar.localize import (extract_regions, get_unique_seeds,
                             unique_seed_string)
from kevlar.tests import data_file


def test_seed_match_set_simple():
    intervals = SeedMatchSet(seedsize=25)
    assert intervals.get_spans('chr1') is None

    intervals.add('chr1', 100)
    intervals.add('chr1', 115)
    assert intervals.get_spans('chr1') == [(100, 140)]

    intervals.add('chr2', 200)
    intervals.add('chr2', 205)
    intervals.add('chr2', 207)
    intervals.add('chr2', 235008)
    intervals.add('chr2', 235075)
    assert intervals.get_spans('chr2') == [(200, 232), (235008, 235100)]


@pytest.mark.parametrize('infile', [
    ('smallseq.fa.gz'),
    ('smallseq.augfasta'),
    ('smallseq.fq'),
])
def test_get_unique_seeds(infile):
    infile = data_file(infile)
    instream = kevlar.parse_augmented_fastx(kevlar.open(infile, 'r'))
    seeds = set([k for k in get_unique_seeds(instream, seedsize=9)])
    testseeds = set(
        ['TTAATTGGC', 'CTTAATTGG', 'TAATTGGCC', 'ATTACCGGT',
         'TTACCGGTA', 'CCTTAATTG', 'GCCTTAATT', 'GGCCTTAAT']
    )
    print(sorted(seeds))
    print(sorted(testseeds))
    assert seeds == testseeds


def test_unique_seed_string():
    infile = data_file('smallseq.augfasta')
    instream = kevlar.parse_augmented_fastx(kevlar.open(infile, 'r'))
    fastastring = unique_seed_string(instream, seedsize=9)
    fastafile = StringIO(fastastring)
    seeds = set([s for d, s in kevlar.seqio.parse_fasta(fastafile)])
    testseeds = set(
        ['TTAATTGGC', 'CTTAATTGG', 'TAATTGGCC', 'ATTACCGGT',
         'TTACCGGTA', 'CCTTAATTG', 'GCCTTAATT', 'GGCCTTAAT']
    )
    print(sorted(seeds))
    print(sorted(testseeds))
    assert seeds == testseeds


def test_extract_regions_basic():
    intervals = SeedMatchSet(seedsize=10)
    intervals.add('bogus-genome-chr2', 10)
    instream = open(data_file('bogus-genome/refr.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=0)]
    assert len(regions) == 1
    assert regions[0] == ('bogus-genome-chr2_10-20', 'GTTACATTAC')


def test_extract_regions_basic_2():
    intervals = SeedMatchSet(seedsize=21)
    intervals.add('simple', 49)
    intervals.add('simple', 52)
    intervals.add('simple', 59)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=5)]
    assert len(regions) == 1
    assert regions[0] == ('simple_44-85',
                          'AATACTATGCCGATTTATTCTTACACAATTAAATTGCTAGT')


def test_extract_regions_large_span():
    intervals = SeedMatchSet(seedsize=21)
    intervals.add('simple', 100)
    intervals.add('simple', 200)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    seqids = [r[0] for r in extract_regions(instream, intervals, maxdiff=50)]
    assert seqids == ['simple_75-146', 'simple_175-246']

    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=50)]
    assert len(regions) == 1
    assert regions[0][0] == 'simple_50-271'


def test_extract_regions_missing_seq():
    intervals = SeedMatchSet(seedsize=21)
    intervals.add('simple', 100)
    intervals.add('simple', 200)
    intervals.add('TheCakeIsALie', 42)
    intervals.add('TheCakeIsALie', 100)
    intervals.add('TheCakeIsALie', 77)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    with pytest.raises(KevlarRefrSeqNotFoundError) as rnf:
        _ = [r for r in extract_regions(instream, intervals)]
    assert 'TheCakeIsALie' in str(rnf)


def test_extract_regions_boundaries():
    intervals = SeedMatchSet(seedsize=31)
    intervals.add('simple', 15)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=20)]
    assert len(regions) == 1
    assert regions[0][0] == 'simple_0-66'

    intervals = SeedMatchSet(seedsize=31)
    intervals.add('simple', 925)
    intervals.add('simple', 955)
    intervals.add('simple', 978)
    instream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    regions = [r for r in extract_regions(instream, intervals, delta=20)]
    assert len(regions) == 1
    assert regions[0][0] == 'simple_905-1000'


def test_main(capsys):
    contig = data_file('localize-contig.fa')
    refr = data_file('localize-refr.fa')
    arglist = ['localize', '--seed-size', '23', '--delta', '50', contig, refr]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.localize.main(args)
    out, err = capsys.readouterr()
    assert '>seq1_10-191' in out


def test_main_no_matches(capsys):
    contig = data_file('localize-contig-bad.fa')
    refr = data_file('localize-refr.fa')
    arglist = ['localize', '--seed-size', '23', contig, refr]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.localize.main(args)
    out, err = capsys.readouterr()
    assert 'WARNING: no reference matches' in err
