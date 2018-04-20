#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import sys
import kevlar
from kevlar.tests import data_file


@pytest.mark.parametrize('greedy', [True, False])
def test_pico_4(greedy, capsys):
    reads = data_file('pico-4.augfastq.gz')
    refr = data_file('human-random-pico.fa.gz')
    arglist = ['alac', '--ksize', '25', reads, refr]
    args = kevlar.cli.parser().parse_args(arglist)
    args.greedy = greedy
    kevlar.alac.main(args)
    out, err = capsys.readouterr()

    vcf = '\t'.join([
        'seq1', '1175768', '.', 'T', 'C', '.', 'PASS',
        'ALTWINDOW=CCCTGCCATTATAGATGCTAGATTCACATCTTCATTTATTTTTACTTTT;'
        'CIGAR=50D192M50D;IKMERS=25;KSW2=179;'
        'REFRWINDOW=CCCTGCCATTATAGATGCTAGATTTACATCTTCATTTATTTTTACTTTT;'
        'CONTIG=ACCTGATTTTGAAGAAGAAAATCAGTTTAAGTCAAAAGGTTACTTTCCTTGTCCTGAACTGG'
        'AGAACTGGGGCCCTGCCATTATAGATGCTAGATTCACATCTTCATTTATTTTTACTTTTTGTCTTGACA'
        'GAGTGGGCGCTGGTTTTTTTAATTATTTTTGGCCAATCAAAAAATACTCTCCTTCGTGGGT'
    ])
    assert vcf.strip() == out.strip()


@pytest.mark.parametrize('cc,pos,ref,alt', [
    (2, 834645, 'A', 'AGTGGGATTACGTAGGAAATCCGCGGGGCTGTGACATATATTTGTTGACAAGCATA'
                     'TATTGTTCCTAGAGGTCGTTGGGTTCGTTACACCCAAGGGGGCGTATAACATGTTA'
                     'CTCAGTTGCGTCGGACCGATTAATAACTCGAATGTAAGGCAGGATATTT'),
    (3, 4072, 'G', 'GCCGAGACGCAGCGTGATACTTAAGATTAAGTTAAGCAACAGCTTAGCGTACGCAATT'
                   'GCGTCTAATTGAGGGGCCGTAGATATAAGCTCCGTGTTCTCAGTTGGTGGGTAACAGA'
                   'ACCCGCAAGCACACCGCTTTCAGTGTGTCACATGCACA'),
    (5, 1175767, 'T', 'C'),
    (6, 185751, 'TCAAACTCTGGCATTATACATAGGGTTCCCG', 'T'),
    (8, 636698, 'C', 'A'),
    (10, 1527138, 'C', 'CTCCTGGTCTGCCACGGTTGACTTGCCTACATAT'),
])
def test_pico_calls(cc, pos, ref, alt):
    reads = data_file('pico-var/cc{:d}.afq.gz'.format(cc))
    readstream = kevlar.parse_augmented_fastx(kevlar.open(reads, 'r'))
    pstream = kevlar.parse_partitioned_reads(readstream)
    refrfile = data_file('human-random-pico.fa.gz')
    caller = kevlar.alac.alac(pstream, refrfile, ksize=25, delta=50)
    calls = [v for v in caller]

    assert len(calls) == 1
    assert calls[0]._pos == pos
    assert calls[0]._refr == ref
    assert calls[0]._alt == alt


def test_pico_partitioned(capsys):
    reads = data_file('pico-partitioned.augfastq.gz')
    refr = data_file('pico-trio-refr.fa.gz')
    arglist = ['alac', '--delta', '50', reads, refr]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.alac.main(args)

    out, err = capsys.readouterr()
    lines = out.strip().split('\n')
    assert len(lines) == 10
    numnocalls = sum([1 for line in lines if '\t.\t.\t.\t.\t' in line])
    assert numnocalls == 2


def test_ikmer_filter_python():
    """
    Smoke test for filtering based in number of supporting ikmers.

    Each partition in the data set has only 2 supporting interesting k-mers.
    The supplied reference file doesn't actually correspond to the reads, so if
    this test passes it's because the filtering worked correctly and the
    `localize` code is never invoked.
    """
    readfile = data_file('min_ikmers_filt.augfastq.gz')
    reads = kevlar.parse_augmented_fastx(kevlar.open(readfile, 'r'))
    parts = kevlar.parse_partitioned_reads(reads)
    refr = data_file('localize-refr.fa')
    calls = list(kevlar.alac.alac(parts, refr, ksize=31, min_ikmers=3))


def test_ikmer_filter_cli():
    reads = data_file('min_ikmers_filt.augfastq.gz')
    refr = data_file('localize-refr.fa')
    arglist = ['alac', '--ksize', '31', '--min-ikmers', '3', reads, refr]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.alac.main(args)


def test_no_reference_match(capsys):
    readfile = data_file('pico-4.augfastq.gz')
    reads = kevlar.parse_augmented_fastx(kevlar.open(readfile, 'r'))
    partitions = kevlar.parse_partitioned_reads(reads)
    refr = data_file('localize-refr.fa')
    baldwin = kevlar.alac.alac(partitions, refr, logstream=sys.stderr)
    calls = list(baldwin)
    out, err = capsys.readouterr()
    assert 'WARNING: no reference matches' in err


@pytest.mark.parametrize('label,position', [
    ('1', 284801),
    ('2', 1660735),
    ('3', 2315888),
    ('4', 2321205),
    ('5', 593252),
])
def test_alac_single_partition(label, position):
    readfile = data_file('fiveparts.augfastq.gz')
    refrfile = data_file('fiveparts-refr.fa.gz')
    readstream = kevlar.parse_augmented_fastx(kevlar.open(readfile, 'r'))
    partstream = kevlar.parse_single_partition(readstream, label)
    calls = list(kevlar.alac.alac(partstream, refrfile))
    assert len(calls) == 1
    assert calls[0].position == position - 1
    assert calls[0].attribute('PART') == label


def test_alac_single_partition_badlabel(capsys):
    readfile = data_file('fiveparts.augfastq.gz')
    refrfile = data_file('fiveparts-refr.fa.gz')
    arglist = ['alac', '--part-id', '6', readfile, refrfile]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.alac.main(args)
    out, err = capsys.readouterr()
    assert out == ''


def test_alac_bigpart():
    readfile = data_file('fiveparts.augfastq.gz')
    refrfile = data_file('fiveparts-refr.fa.gz')
    readstream = kevlar.parse_augmented_fastx(kevlar.open(readfile, 'r'))
    partstream = kevlar.parse_partitioned_reads(readstream)
    calls = list(kevlar.alac.alac(partstream, refrfile, bigpart=20))
    assert len(calls) == 3


@pytest.mark.parametrize('cc,numcalls,numdist,numinf', [
    ('26849', 4, 2, 0),
    ('138713', 14, 14, 10),
])
def test_alac_inf_mate_dist(cc, numcalls, numdist, numinf):
    readfile = data_file('inf-mate-dist/cc{}.augfastq.gz'.format(cc))
    refrfile = data_file('inf-mate-dist/cc{}.genome.fa.gz'.format(cc))
    readstream = kevlar.parse_augmented_fastx(kevlar.open(readfile, 'r'))
    partstream = kevlar.parse_partitioned_reads(readstream)
    caller = kevlar.alac.alac(partstream, refrfile, ksize=31, delta=50,
                              seedsize=51, fallback=True)
    calls = list(caller)
    withdist = [c for c in calls if c.attribute('MD') is not None]
    infdist = [c for c in calls if c.attribute('MD') == 'inf']
    print('DEBUG total withdist infdist', len(calls), len(withdist),
          len(infdist), file=sys.stderr)
    assert len(calls) == numcalls
    assert len(withdist) == numdist
    assert len(infdist) == numinf


def test_alac_no_mates():
    readfile = data_file('inf-mate-dist/cc26849-nomates.augfastq.gz')
    refrfile = data_file('inf-mate-dist/cc26849.genome.fa.gz')
    readstream = kevlar.parse_augmented_fastx(kevlar.open(readfile, 'r'))
    partstream = kevlar.parse_partitioned_reads(readstream)
    caller = kevlar.alac.alac(partstream, refrfile, ksize=31, delta=50,
                              seedsize=51, fallback=True)
    calls = list(caller)
    print(*[c.vcf for c in calls], sep='\n', file=sys.stderr)
    assert len(calls) == 4
    filtcalls = [c for c in calls if c.filterstr == 'PASS']
    assert len(filtcalls) == 3


@pytest.mark.parametrize('vcfposition,X,cigar', [
    (40692, 10000, '32713D96M6I91M15142D'),
    (40692, 1000, '50D96M6I91M50D'),
    (40692, 0, '32713D96M6I91M140025D'),
    (40692, None, '50D96M6I91M50D'),
])
def test_alac_maxdiff(vcfposition, X, cigar):
    pstream = kevlar.parse_partitioned_reads(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file('maxdiff-reads.augfastq.gz'), 'r')
        )
    )
    refrfile = data_file('maxdiff-refr.fa.gz')
    caller = kevlar.alac.alac(
        pstream, refrfile, ksize=31, delta=50, seedsize=51, maxdiff=X
    )
    calls = list(caller)
    assert len(calls) == 1
    assert calls[0].cigar == cigar
    assert calls[0].position == vcfposition - 1
