#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019 Battelle National Biodefense Institute
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from intervaltree import Interval
import kevlar
from kevlar.tests import data_file
import pytest


def test_load_variant_mask():
    bedfile = data_file('fiveparts-ignore-single.bed')
    with kevlar.open(bedfile, 'r') as bedstream:
        index = kevlar.varfilter.load_variant_mask(bedstream)
    assert len(index) == 1
    assert list(index.trees.keys()) == ['chr17']
    assert index.query('chr1', 1, 10000000) == set()
    assert index.query('chr17', 1, 10000000) == set()
    itvl = Interval(36385017, 36385018, 'chr17:36385017-36385018')
    assert index.query('chr17', 36385017) == set([itvl])


def test_load_variant_mask_multi_chrom():
    bedfile = data_file('intervals.bed')
    with kevlar.open(bedfile, 'r') as bedstream:
        index = kevlar.varfilter.load_variant_mask(bedstream)
    assert len(index) == 5
    assert list(index.trees.keys()) == ['1', '2', '3', '4']
    assert index.query('5', 1, 10000000) == set()
    itvls = [
        Interval(40, 400, '4:40-400'),
        Interval(44, 444, '4:44-444'),
    ]
    assert index.query('4', 300) == set(itvls)


def test_varfilter_single():
    bedfile = data_file('fiveparts-ignore-single.bed')
    with kevlar.open(bedfile, 'r') as bedstream:
        index = kevlar.varfilter.load_variant_mask(bedstream)
    vcffile = data_file('five-snvs-with-likelihood.vcf')
    with kevlar.open(vcffile, 'r') as vcfstream:
        reader = kevlar.vcf.VCFReader(vcfstream)
        varcalls = list(kevlar.varfilter.varfilter(reader, index))
    assert len(varcalls) == 5
    filtered = [vc for vc in varcalls if vc.filterstr != 'PASS']
    assert len(filtered) == 1
    assert filtered[0].position == 36385017


def test_varfilter_main(capsys):
    arglist = [
        'varfilter', data_file('fiveparts-ignore.bed'),
        data_file('five-snvs-with-likelihood.vcf')
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.varfilter.main(args)

    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    calls = [line for line in outlines if not line.startswith('#')]
    assert len(calls) == 5
    filtered = [c for c in calls if '\tUserFilter\t' in c]
    assert len(filtered) == 2
    assert [c.split('\t')[1] for c in filtered] == ['36385018', '3547691']
