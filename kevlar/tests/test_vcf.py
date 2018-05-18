#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import pytest
import kevlar
from kevlar.tests import data_file
from kevlar.vcf import Variant, FormattedList
from kevlar.vcf import VariantFilter as vf


@pytest.fixture
def yrb_writer():
    writer = kevlar.vcf.VCFWriter(sys.stdout, source='py.test')
    writer.register_sample('NA19238')
    writer.register_sample('NA19239')
    writer.register_sample('NA19240')
    writer.describe_format('GT', 'String', '1', 'Genotype')
    return writer


def test_snv_obj():
    snv = Variant('scaffold42', 10773, 'A', 'G')
    assert str(snv) == 'scaffold42:10773:A->G'
    vcfvalues = ['scaffold42', '10774', '.', 'A', 'G', '.', 'PASS', '.']
    assert snv.vcf == '\t'.join(vcfvalues)
    assert snv.cigar is None

    snv2 = Variant('chr5', 500, 'T', 'G', CIGAR='10D200M10D')
    assert snv2.cigar == '10D200M10D'
    assert snv2.window is None


def test_indel_obj():
    """
    Test indel objects

    The coordinate used to construct the object is 0-based, but includes the
    nucleotide shared by the reference and alternate alleles. The str() output
    coordinate is increased by 1 to account for this nucleotide, while the VCF
    output is increased by 1 to transform to a 1-based system where the shared
    nucleotide is the point of reference.
    """
    indel1 = Variant('chr3', 8998622, 'GATTACA', 'G')
    assert str(indel1) == 'chr3:8998623:6D'
    vcfvalues = ['chr3', '8998623', '.', 'GATTACA', 'G', '.', 'PASS', '.']
    assert indel1.vcf == '\t'.join(vcfvalues)

    indel2 = Variant('chr6', 75522411, 'G', 'GATTACA')
    assert str(indel2) == 'chr6:75522412:I->ATTACA'
    vcfvalues = ['chr6', '75522412', '.', 'G', 'GATTACA', '.', 'PASS', '.']
    assert indel2.vcf == '\t'.join(vcfvalues)


def test_filter_field():
    v = Variant('scaffold1', 12345, '.', '.')
    assert v.filterstr == '.'
    v.filter(vf.InscrutableCigar)
    assert v.filterstr == 'InscrutableCigar'

    v = Variant('chr1', 55555, '.', '.')
    v.filter(vf.PerfectMatch)
    assert v.filterstr == 'PerfectMatch'

    v = Variant('1', 809768, 'C', 'CAT')
    assert v.filterstr == 'PASS'
    v.filter(vf.PassengerVariant)
    assert v.filterstr == 'PassengerVariant'
    v.filter(vf.MateFail)
    assert v.filterstr == 'MateFail;PassengerVariant'

    v = Variant('one', 112358, 'T', 'A')
    v.filter('SNPyMcSNPface')
    v.filter(6.022e23)
    v.filter(dict(chicken='waffles', biscuits='gravy'))
    v.filterstr == 'PASS'  # These "filters" shouldn't actually do anything


def test_info():
    """Test handling of "info" field attributes.

    This tests the mechanics of the .annotate() and .attribute() API, and the
    FormattedList class underpinning it.
    """
    values = FormattedList()
    assert str(values) == '.'
    values.append(42)
    assert str(values) == '42'
    values.append(1776)
    assert str(values) == '42,1776'
    values.append('B0gU$')
    with pytest.raises(kevlar.vcf.KevlarMixedDataTypeError):
        str(values)

    v = Variant('1', 12345, 'G', 'C')
    assert v.attribute('VW') is None

    v.annotate('VW', 'GATTACA')
    assert v.attribute('VW') == 'GATTACA'
    assert v.attribute('VW', pair=True) == 'VW=GATTACA'

    v.annotate('VW', 'ATGCCCTAG')
    assert v.info['VW'] == ['GATTACA', 'ATGCCCTAG']
    assert v.attribute('VW') == ['GATTACA', 'ATGCCCTAG']
    assert v.attribute('VW', string=True) == 'GATTACA,ATGCCCTAG'
    assert v.attribute('VW', pair=True) == 'VW=GATTACA,ATGCCCTAG'

    v.annotate('VW', 'AAAAAAAAA')
    assert v.attribute('VW') == ['GATTACA', 'ATGCCCTAG', 'AAAAAAAAA']
    assert v.attribute('VW', pair=True) == 'VW=GATTACA,ATGCCCTAG,AAAAAAAAA'

    v.annotate('DROPPED', 3)
    assert v.attribute('DROPPED') == 3
    assert v.attribute('DROPPED', string=True) == '3'

    v.annotate('DROPPED', 31)
    assert v.attribute('DROPPED') == [3, 31]
    assert v.attribute('DROPPED', string=True) == '3,31'
    assert v.attribute('DROPPED', pair=True) == 'DROPPED=3,31'

    v.annotate('MATEDIST', 432.1234)
    v.annotate('MATEDIST', 8765.4321)
    assert v.attribute('MATEDIST', string=True) == '432.123,8765.432'

    v.annotate('LLIH', -436.0111857750478)
    assert v.attribute('LLIH', pair=True) == 'LLIH=-436.011'


def test_format():
    v = Variant('1', 12345, 'G', 'C')
    v.format('NA19238', 'GT', '0/0')
    assert v.format('NA19238', 'GT') == '0/0'
    assert v.format('NA19238', 'XYZ') is None
    assert v.format('NA19239', 'GT') is None


def test_writer(yrb_writer, capsys):
    yrb_writer = kevlar.vcf.VCFWriter(sys.stdout, source='py.test')
    yrb_writer.register_sample('NA19238')
    yrb_writer.register_sample('NA19239')
    yrb_writer.register_sample('NA19240')
    yrb_writer.describe_format('GT', 'String', '1', 'Genotype')
    yrb_writer.write_header()

    v = Variant('1', 12345, 'G', 'C')
    v.annotate('PART', '42')
    v.annotate('CONTIG', 'A' * 100)
    v.format('NA19238', 'GT', '0/0')
    v.format('NA19239', 'GT', '0/0')
    v.format('NA19240', 'GT', '0/1')
    v.format('NA19238', 'ALTABUND', '12,9,8')
    v.format('NA19239', 'ALTABUND', '0,0,0')
    v.format('NA19240', 'ALTABUND', '0,0,0')
    yrb_writer.write(v)

    out, err = capsys.readouterr()
    print(out)

    outlines = out.strip().split('\n')
    fmtlines = [l for l in outlines if l.startswith('##FORMAT')]
    assert len(fmtlines) == 2
    gtfmt = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    assert gtfmt in fmtlines

    varlines = [l for l in outlines if not l.startswith('#')]
    assert len(varlines) == 1
    values = varlines[0].split('\t')
    assert len(values) == 12
    assert values[8:12] == [
        'ALTABUND:GT', '12,9,8:0/0', '0,0,0:0/0', '0,0,0:0/1'
    ]


def test_writer_bad_fmt(yrb_writer):
    v = Variant('1', 12345, 'G', 'C')
    v.annotate('PART', '42')
    v.annotate('CONTIG', 'A' * 100)
    v.format('NA19238', 'GT', '0/0')
    v.format('NA19240', 'GT', '0/1')
    v.format('NA19239', 'ALTABUND', '0,0,0')
    v.format('NA19240', 'ALTABUND', '0,0,0')
    with pytest.raises(kevlar.vcf.VariantAnnotationError) as vae:
        yrb_writer.write(v)
    assert 'samples not annotated with the same FORMAT fields' in str(vae)


def test_reader():
    instream = kevlar.open(data_file('five-snvs-with-likelihood.vcf'), 'r')
    reader = kevlar.vcf.VCFReader(instream)
    calls = list(reader)
    print(calls)
    assert len(calls) == 5
    assert calls[1].attribute('PART') == '54'
    assert calls[3].format('Kid', 'ALTABUND') == (
        '21,20,20,19,17,19,20,19,18,17,17,17,17,17,17,17,18,19,19,19,18,18,18,'
        '17,19,18,17,17,17,15,15'
    )


@pytest.mark.parametrize('filename,errormsg', [
    ('five-snvs-fmt-mismatch.vcf', 'sample number mismatch'),
    ('five-snvs-fmtstr-mismatch.vcf', 'format data mismatch'),
])
def test_reader_format_mismatch(filename, errormsg):
    instream = kevlar.open(data_file(filename), 'r')
    reader = kevlar.vcf.VCFReader(instream)
    with pytest.raises(kevlar.vcf.VariantAnnotationError) as vae:
        calls = list(reader)
    assert errormsg in str(vae)


def test_vcf_roundtrip(capsys):
    instream = kevlar.open(data_file('five-snvs-with-likelihood.vcf'), 'r')
    reader = kevlar.vcf.VCFReader(instream)

    writer = kevlar.vcf.VCFWriter(
        sys.stdout, source=None,
        refr='GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'
    )
    writer.register_sample('Kid')
    writer.register_sample('Mom')
    writer.register_sample('Dad')
    writer.describe_format('GT', 'String', '1', 'Genotype')
    writer.write_header(skipdate=True)
    calls = list()
    for call in reader:
        calls.append(call)
        writer.write(call)

    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    reader2 = kevlar.vcf.VCFReader(outlines)
    calls2 = list(reader2)
    assert len(calls) == len(calls2)
    assert [c.position for c in calls] == [c.position for c in calls2]
    assert [str(c) for c in calls] == [str(c) for c in calls2]
    assert [c.window for c in calls] == [c.window for c in calls2]
