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


@pytest.mark.parametrize('seq,pos,offset,refr,alt,refrwindow,altwindow', [
    ('AACTAGCCTGCGGTCTGTGTTTCCCGACTTCTGAGTCATGGGGTTTCAATGCCTAT',
     14, 2, 'C', 'T', 'CCTGCGGTCTGTGTTTC', 'CCTGCGGTTTGTGTTTC'),
    ('TTGAGATCGCGACGCTACTCTGAGCTCGGAGGAGCGGCATAAACGCGCCACCACCC',
     26, 1, 'C', 'G', 'TCTGAGCTCGGAGGAGC', 'TCTGAGCTGGGAGGAGC'),
    ('CCTTGGTGCCACGATCCGGCTATGGCGGAAGGGCACACCTAACCGCAACATTTGCC',
     2, 2, 'T', 'C', 'CCTTGGTGCCA', 'CCCTGGTGCCA'),
    ('GGGTCCCAAGAGTCTGATTTCTAGCTTTTTATTTACACCCCGGTAGCAGGATCAGA',
     33, 3, 'T', 'G', 'TTTTTATTTACACCCCG', 'TTTTTATTGACACCCCG'),
])
def test_snv(seq, pos, offset, refr, alt, refrwindow, altwindow):
    testrefr, testalt, testrw, testaw = kevlar.gentrio.mutate_snv(
        seq, pos, offset, ksize=9
    )
    
    print('REFR', refr, testrefr, refrwindow, testrw, file=sys.stderr)
    print('ALT', alt, testalt, altwindow, testaw, file=sys.stderr)

    assert refr == testrefr
    assert alt == testalt
    assert refrwindow == testrw
    assert altwindow == testaw


@pytest.mark.parametrize('seq,pos,length,duplpos,refr,alt,rwindow,awindow', [
    ('AACTAGCCTGCGGTCTGTGTTTCCCGACTTCTGAGTCATGGGGTTTCAATGCCTAT',
     11, 5, 33, 'C', 'CAGTCA', 'CTGCGGTC', 'CTGCAGTCAGGTC'),
    ('TTGAGATCGCGACGCTACTCTGAGCTCGGAGGAGCGGCATAAACGCGCCACCACCC',
     47, 11, 32, 'G', 'GAGCGGCATAAA', 'CGCGCCAC', 'CGCGAGCGGCATAAACCAC'),
    ('CCTTGGTGCCACGATCCGGCTATGGCGGAAGGGCACACCTAACCGCAACATTTGCC',
     52, 3, 39, 'T', 'TTAA', 'CATTTGCC', 'CATTTAATGCC'),
    ('GGGTCCCAAGAGTCTGATTTCTAGCTTTTTATTTACACCCCGGTAGCAGGATCAGA',
     9, 9, 29, 'A', 'ATATTTACAC', 'CCAAGAGT', 'CCAATATTTACACGAGT'),
])
def test_insertion(seq, pos, length, duplpos, refr, alt, rwindow, awindow):
    testrefr, testalt, testrw, testaw = kevlar.gentrio.mutate_insertion(
        seq, pos, length, duplpos, ksize=5
    )
    
    print('REFR', refr, testrefr, rwindow, testrw, file=sys.stderr)
    print('ALT', alt, testalt, awindow, testaw, file=sys.stderr)

    assert refr == testrefr
    assert alt == testalt
    assert rwindow == testrw
    assert awindow == testaw