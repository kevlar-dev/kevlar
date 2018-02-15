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
from kevlar.localize import Localizer, KevlarRefrSeqNotFoundError
from kevlar.localize import get_unique_seeds, unique_seed_string
from kevlar.reference import ReferenceCutout
from kevlar.tests import data_file


def test_localizer_simple():
    intervals = Localizer(seedsize=25)
    assert list(intervals.get_cutouts()) == []

    intervals.add_seed_match('chr1', 100)
    intervals.add_seed_match('chr1', 115)
    intervals.add_seed_match('chr2', 200)
    intervals.add_seed_match('chr2', 205)
    intervals.add_seed_match('chr2', 207)
    intervals.add_seed_match('chr2', 235008)
    intervals.add_seed_match('chr2', 235075)
    testint = [c.interval for c in intervals.get_cutouts()]
    print('DEBUG', testint)
    assert testint == [
        ('chr1', 100, 140),
        ('chr2', 200, 232),
        ('chr2', 235008, 235100)
    ]


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


def test_get_cutouts_basic():
    intervals = Localizer(seedsize=10)
    intervals.add_seed_match('bogus-genome-chr2', 10)
    refrstream = open(data_file('bogus-genome/refr.fa'), 'r')
    seqs = kevlar.seqio.parse_seq_dict(refrstream)
    cutouts = list(intervals.get_cutouts(refrseqs=seqs))
    assert len(cutouts) == 1
    assert cutouts[0].defline == 'bogus-genome-chr2_10-20'
    assert cutouts[0].sequence == 'GTTACATTAC'


def test_get_cutouts_basic_2():
    intervals = Localizer(seedsize=21, delta=5)
    intervals.add_seed_match('simple', 49)
    intervals.add_seed_match('simple', 52)
    intervals.add_seed_match('simple', 59)
    refrstream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    seqs = kevlar.seqio.parse_seq_dict(refrstream)
    cutouts = list(intervals.get_cutouts(refrseqs=seqs))
    assert len(cutouts) == 1
    assert cutouts[0].defline == 'simple_44-85'
    assert cutouts[0].sequence == 'AATACTATGCCGATTTATTCTTACACAATTAAATTGCTAGT'


def test_get_cutouts_large_span():
    refrstream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    seqs = kevlar.seqio.parse_seq_dict(refrstream)

    intervals = Localizer(seedsize=21, delta=25)
    intervals.add_seed_match('simple', 100)
    intervals.add_seed_match('simple', 200)
    cutouts = intervals.get_cutouts(refrseqs=seqs, clusterdist=50)
    deflines = [c.defline for c in cutouts]
    assert deflines == ['simple_75-146', 'simple_175-246']

    intervals = Localizer(seedsize=21, delta=50)
    intervals.add_seed_match('simple', 100)
    intervals.add_seed_match('simple', 200)
    cutouts = intervals.get_cutouts(refrseqs=seqs, clusterdist=100)
    deflines = [c.defline for c in cutouts]
    assert deflines == ['simple_50-271']


def test_get_cutouts_missing_seq():
    intervals = Localizer(seedsize=21)
    intervals.add_seed_match('simple', 100)
    intervals.add_seed_match('simple', 200)
    intervals.add_seed_match('TheCakeIsALie', 42)
    intervals.add_seed_match('TheCakeIsALie', 100)
    intervals.add_seed_match('TheCakeIsALie', 77)
    refrstream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    seqs = kevlar.seqio.parse_seq_dict(refrstream)
    with pytest.raises(KevlarRefrSeqNotFoundError) as rnf:
        list(intervals.get_cutouts(refrseqs=seqs))
    assert 'TheCakeIsALie' in str(rnf)


def test_extract_regions_boundaries():
    refrstream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    seqs = kevlar.seqio.parse_seq_dict(refrstream)

    intervals = Localizer(seedsize=31, delta=20)
    intervals.add_seed_match('simple', 15)
    cutouts = list(intervals.get_cutouts(refrseqs=seqs))
    assert len(cutouts) == 1
    assert cutouts[0].defline == 'simple_0-66'

    intervals = Localizer(seedsize=31, delta=20)
    intervals.add_seed_match('simple', 925)
    intervals.add_seed_match('simple', 955)
    intervals.add_seed_match('simple', 978)
    cutouts = list(intervals.get_cutouts(refrseqs=seqs))
    assert len(cutouts) == 1
    assert cutouts[0].defline == 'simple_905-1000'


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
