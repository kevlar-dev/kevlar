#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import kevlar
from kevlar.varmap import VariantMapping
from kevlar.tests import data_file
import pytest
import screed


def test_call_ssc_1bpdel():
    """Test 1bp deletion"""
    qfile = data_file('ssc218.contig.augfasta')
    tfile = data_file('ssc218.gdna.fa')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    query = [record for record in qinstream][0]
    tinstream = kevlar.reference.load_refr_cutouts(kevlar.open(tfile, 'r'))
    target = [record for record in tinstream][0]
    aln = VariantMapping(query, target, 1e6, '50D132M1D125M50D')
    variants = list(aln.call_variants(31))

    assert len(variants) == 1
    assert str(variants[0]) == '6:23230160:1D'


def test_call_ssc_two_proximal_snvs():
    """Test two proximal SNVs

    Currently this serves as a negative control for calling isolated SNVs, but
    distinguishing which (if any) of a set of proximal SNVs is novel will be
    supported soon, and this test will need to be updated.
    """
    qfile = data_file('ssc107.contig.augfasta.gz')
    tfile = data_file('ssc107.gdna.fa.gz')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    query = [record for record in qinstream][0]
    tinstream = kevlar.reference.load_refr_cutouts(kevlar.open(tfile, 'r'))
    target = [record for record in tinstream][0]

    aln = VariantMapping(query, target, 1e6, '25D263M25D')
    variants = list(aln.call_variants(31))
    assert len(variants) == 2


@pytest.mark.parametrize('prefix,cigar,refrwindow,altwindow', [
    ('phony-snv-01', '25D98M25D',
        'GGGGGTGTCTGCGACCACAGCTGAACATGACGAAACGGGTG',
        'GGGGGTGTCTGCGACCACAGGTGAACATGACGAAACGGGTG'),
    ('phony-snv-02', '24D99M25D',
        'ATTCGTATTACCCCTGGGATTTGGGAGCTGGTCTATATAGG',
        'ATTCGTATTACCCCTGGGATATGGGAGCTGGTCTATATAGG'),
    ('phony-deletion-01', '25D28M8D49M25D',
        'GGCTCAAGACTAAAAAGACTGAGACTCGTTTTTGGTGACAAGCAGGGC',
        'GGCTCAAGACTAAAAAGACTTTTTTGGTGACAAGCAGGGC'),
    ('phony-deletion-02', '40D29M3D36M40D',
        'CATCATCTCGTAGGTTTGTCTAGTGCAAACAGAGTCCCCCTGC',
        'CATCATCTCGTAGGTTTGTCTGCAAACAGAGTCCCCCTGC'),
    ('phony-insertion-01', '10D34M7I49M10D1M',
        'CATCTGTTTTTCTCGAACTCGTATATTATCTATAAATTCC',
        'CATCTGTTTTTCTCGAACTCGATTACAGTATATTATCTATAAATTCC'),
    ('phony-insertion-02', '10D33M27I95M10D',
        'GCCAGGAAGTTTACGATAAGGTGTTGCCATTCGAAATGAC',
        'GCCAGGAAGTTTACGATAAGTATATATATATATATATATATATATATGTGTTGCCATTCGAAATGAC'),
])
def test_variant_window(prefix, cigar, refrwindow, altwindow):
    qfile = data_file(prefix + '.contig.fa')
    tfile = data_file(prefix + '.gdna.fa')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    query = [record for record in qinstream][0]
    tinstream = kevlar.reference.load_refr_cutouts(kevlar.open(tfile, 'r'))
    target = [record for record in tinstream][0]

    aln = VariantMapping(query, target, 1e6, cigar)
    variants = list(aln.call_variants(21))
    assert len(variants) == 1
    assert variants[0].window == altwindow
    assert variants[0].refrwindow == refrwindow


def test_nocall():
    # Intentionally mismatched
    qfile = data_file('phony-deletion-01.contig.fa')
    tfile = data_file('phony-insertion-01.gdna.fa')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    query = [record for record in qinstream][0]
    tinstream = kevlar.reference.load_refr_cutouts(kevlar.open(tfile, 'r'))
    target = [record for record in tinstream][0]

    aln = VariantMapping(query, target, 1e6, '25D5M22I5M46D8M13D2M35I')
    assert aln.offset is None
    assert aln.targetshort is None
    assert aln.match is None
    assert aln.leftflank is None
    assert aln.indel is None
    assert aln.indeltype is None
    assert aln.rightflank is None

    variants = list(aln.call_variants(21))
    assert len(variants) == 1
    assert variants[0].vcf == (
        'yourchr\t801\t.\t.\t.\t.\tInscrutableCigar\t'
        'CIGAR=25D5M22I5M46D8M13D2M35I;KSW2=1000000.0;CONTIG=AACTGGTGGGCTCAAGA'
        'CTAAAAAGACTTTTTTGGTGACAAGCAGGGCGGCCTGCCCTTCCTGTAGTGCAAGAAAAT'
    )


def test_variant_mapping():
    contig = screed.Record(
        name='contig1',
        sequence='CCTGAGCCCTCTCAAGTCGGGTCCTGGCCCGGTCTGCCCATGAGGCTGGGCCTGAGCCCC'
    )
    cutout = kevlar.reference.ReferenceCutout(
        defline='chr1_10000-10060',
        sequence='CCTGAGCCCTCTCAAGTCGGGTCCTGGCCCAGTCTGCCCATGAGGCTGGGCCTGAGCCCC'
    )
    mapping = VariantMapping(contig, cutout, score=1e6, cigar='60M')

    assert mapping.seqid == 'chr1'
    assert mapping.interval == ('chr1', 10000, 10060)


@pytest.mark.parametrize('query,target,dist,n,trimcount', [
    ('phony-snv-01b.contig.fa', 'phony-snv-01.gdna.fa', 5, 1, 1),
    ('phony-snv-02b.contig.fa', 'phony-snv-02.gdna.fa', 5, 1, 1),
    ('phony-snv-01b.contig.fa', 'phony-snv-01.gdna.fa', 2, 2, 0),
    ('phony-snv-02b.contig.fa', 'phony-snv-02.gdna.fa', None, 2, 0),
])
def test_call_near_end(query, target, dist, n, trimcount):
    contig = next(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file(query), 'r')
        )
    )
    cutout = next(
        kevlar.reference.load_refr_cutouts(
            kevlar.open(data_file(target), 'r')
        )
    )
    aln = VariantMapping(contig, cutout)
    calls = list(aln.call_variants(31, mindist=dist))
    assert len(calls) == n
    assert aln.trimmed == trimcount


@pytest.mark.parametrize('query,target,vw,rw', [
    (
        'trunc-snv.contig.fa', 'trunc-snv.gdna.fa',
        'TAGCATACAGGTAGTCAGGGGGTGTCTGCGACCACAGCTGAA',
        'TAGCATACAGGAAGTCAGGGGGTGTCTGCGACCACAGCTGAA'
    ),
    (
        'trunc-snv-funky.contig.fa', 'trunc-snv-funky.gdna.fa',
        'TGTGTCTGAGAGGGTGTTGCCAAAGGAGATTAACATTTG',
        'TGTGTCTGTGAGGGTGTTGCCAAAGGAGATTAACATTTG'
    ),
    (
        'trunc-indel-funky.contig.fa', 'trunc-snv-funky.gdna.fa',
        'TGTGTCTGTGAGTATATAGGTGTTGCCAAAGGAGATTAACATTTGAGT',
        'TGTGTCTGTGAGGGTGTTGCCAAAGGAGATTAACATTTGAGT'
    ),
])
def test_call_truncated_windows(query, target, vw, rw):
    contig = next(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file(query), 'r')
        )
    )
    cutout = next(
        kevlar.reference.load_refr_cutouts(
            kevlar.open(data_file(target), 'r')
        )
    )
    aln = VariantMapping(contig, cutout)
    if aln.vartype == 'snv':
        assert aln.leftflank is None
        assert aln.indeltype is None
        assert aln.indel is None
        assert aln.rightflank is None

    calls = list(aln.call_variants(31))
    assert len(calls) == 1
    print('VW:', calls[0].window, file=sys.stderr)
    print('RW:', calls[0].refrwindow, file=sys.stderr)
    assert calls[0].window == vw
    assert calls[0].refrwindow == rw


def test_call_indel_snv():
    contig = next(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file('indel-snv.contig.augfasta'), 'r')
        )
    )
    cutout = next(
        kevlar.reference.load_refr_cutouts(
            kevlar.open(data_file('indel-snv.gdna.fa'), 'r')
        )
    )
    aln = VariantMapping(contig, cutout)
    calls = list(aln.call_variants(31))
    assert len(calls) == 2

    assert calls[0]._refr == 'CA'
    assert calls[0]._alt == 'C'
    assert calls[0]._pos == 501 - 1

    assert calls[1]._refr == 'C'
    assert calls[1]._alt == 'A'
    assert calls[1]._pos == 474 - 1

    calls = list(aln.call_variants(31, mindist=None))
    assert len(calls) == 2


def test_call_num_interesting_kmers():
    contig = next(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file('iktest.contig.fa'), 'r')
        )
    )
    cutout = next(
        kevlar.reference.load_refr_cutouts(
            kevlar.open(data_file('iktest.gdna.fa'), 'r')
        )
    )
    aln = VariantMapping(contig, cutout)
    calls = list(aln.call_variants(29))
    assert len(calls) == 1
    assert calls[0].attribute('IKMERS') == '1'


def test_passenger_screen():
    contig = next(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file('wasp-pass.contig.augfasta'), 'r')
        )
    )
    cutout = next(
        kevlar.reference.load_refr_cutouts(
            kevlar.open(data_file('wasp.gdna.fa'), 'r')
        )
    )
    aln = VariantMapping(contig, cutout)
    calls = list(aln.call_variants(29))
    assert len(calls) == 2
    assert calls[0].filterstr == 'PASS'
    assert calls[1].filterstr == 'PassengerVariant'


@pytest.mark.parametrize('query,target,refr,alt', [
    ('nomargin-snv-contigs.augfasta', 'nomargin-gdna.fa', 'A', 'G'),
    ('nomargin-indel-contigs.augfasta', 'nomargin-gdna.fa', 'AAGT', 'A'),
    ('nomargin-r-snv-contigs.augfasta', 'nomargin-r-gdna.fa', 'A', 'G'),
    ('nomargin-r-indel-contigs.augfasta', 'nomargin-r-gdna.fa', 'C', 'CTAT'),
])
def test_no_margin(query, target, refr, alt):
    contig = next(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file(query), 'r')
        )
    )
    cutout = next(
        kevlar.reference.load_refr_cutouts(
            kevlar.open(data_file(target), 'r')
        )
    )
    aln = VariantMapping(contig, cutout)
    calls = list(aln.call_variants(31))
    assert len(calls) == 1
    assert calls[0].filterstr == 'PASS'
    assert calls[0]._refr == refr
    assert calls[0]._alt == alt


def test_varmap_str():
    contig = next(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file('wasp-pass.contig.augfasta'), 'r')
        )
    )
    cutout = next(
        kevlar.reference.load_refr_cutouts(
            kevlar.open(data_file('wasp.gdna.fa'), 'r')
        )
    )
    aln = VariantMapping(contig, cutout)
    alignstr = kevlar.open(data_file('wasp-align.txt'), 'r').read().strip()

    print(str(aln), file=sys.stderr)
    print(alignstr, file=sys.stderr)

    assert str(aln) == alignstr


def test_drop_numerous_mismatches():
    contig = next(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file('drop-polysnp-contig.augfasta'), 'r')
        )
    )
    cutout = next(
        kevlar.reference.load_refr_cutouts(
            kevlar.open(data_file('drop-polysnp-gdna.fa'), 'r')
        )
    )
    aln = VariantMapping(contig, cutout)
    calls = list(aln.call_variants(21))
    for c in calls:
        print(c.vcf)
    assert len(calls) == 1
    assert calls[0].filterstr == 'NumerousMismatches'
    assert calls[0]._refr == '.'
    assert calls[0]._alt == '.'
