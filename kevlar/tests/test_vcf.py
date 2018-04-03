#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
from kevlar.vcf import Variant, VariantIndel, VariantSNV
from kevlar.vcf import VariantFilter as vf


def test_snv_obj():
    snv = VariantSNV('scaffold42', 10773, 'A', 'G')
    assert str(snv) == 'scaffold42:10773:A->G'
    vcfvalues = ['scaffold42', '10774', '.', 'A', 'G', '.', 'PASS', '.']
    assert snv.vcf == '\t'.join(vcfvalues)
    assert snv.cigar is None

    snv2 = VariantSNV('chr5', 500, 'T', 'G', CG='10D200M10D')
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
    indel1 = VariantIndel('chr3', 8998622, 'GATTACA', 'G')
    assert str(indel1) == 'chr3:8998623:6D'
    vcfvalues = ['chr3', '8998623', '.', 'GATTACA', 'G', '.', 'PASS', '.']
    assert indel1.vcf == '\t'.join(vcfvalues)

    indel2 = VariantIndel('chr6', 75522411, 'G', 'GATTACA')
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
