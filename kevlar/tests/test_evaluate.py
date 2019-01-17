#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 Battelle National Biodefense Institute
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import kevlar
from kevlar.tests import data_file


def test_compact():
    bedstream = kevlar.open(data_file('compact-test-refr.bed.gz'), 'r')
    index = kevlar.evaluate.populate_index_from_bed(bedstream)

    vcfstream = kevlar.open(data_file('compact-test-pred.vcf.gz'), 'r')
    reader = kevlar.vcf.VCFReader(vcfstream)
    compactor = kevlar.evaluate.compact(reader, index, delta=10)
    calls = list(compactor)
    assert len(calls) == 33
