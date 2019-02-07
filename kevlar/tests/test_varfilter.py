#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019 Battelle National Biodefense Institute
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
from kevlar.tests import data_file


def test_load_predictions():
    with kevlar.open(data_file('five-snvs-with-likelihood.vcf'), 'r') as vcf:
        vcfreader = kevlar.vcf.VCFReader(vcf)
        index = kevlar.varfilter.load_predictions(vcfreader)
    assert len(index) == 5
    assert list(index.trees.keys()) == ['chr17']
    assert index.query('chr1', 1, 1000000) == set()
    assert index.query('chr17', 1, 1000000) == set()
    result = [i.data.region for i in index.query('chr17', 36385017)]
    assert result == [('chr17', 36385017, 36385018)]


def test_load_predictions_multi_chrom():
    with kevlar.open(data_file('case-low-abund/calls.vcf.gz'), 'r') as vcf:
        vcfreader = kevlar.vcf.VCFReader(vcf)
        index = kevlar.varfilter.load_predictions(vcfreader)
    assert len(index) == 5
    assert set(index.trees.keys()) == set(['1', '9', '14'])
    assert index.query('chr1', 1, 1000000) == set()
    assert index.query('1', 1, 1000000) == set()
    result = [i.data.region for i in index.query('1', 91850000, 91860000)]
    assert set(result) == set([
        ('1', 91853096, 91853097),
        ('1', 91853110, 91853111),
    ])
    result = [i.data.region for i in index.query('14', 82461000, 82462000)]
    assert result == [('14', 82461856, 82461857)]


def test_varfilter_single():
    bedstream = kevlar.parse_bed(
        kevlar.open(data_file('fiveparts-ignore-single.bed'), 'r')
    )
    vcffile = data_file('five-snvs-with-likelihood.vcf')
    with kevlar.open(vcffile, 'r') as vcfstream:
        reader = kevlar.vcf.VCFReader(vcfstream)
        varcalls = list(kevlar.varfilter.varfilter(reader, bedstream))
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
    positions = [c.split('\t')[1] for c in filtered]
    assert sorted(positions) == sorted(['36385018', '3547691'])
