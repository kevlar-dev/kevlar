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


@pytest.mark.parametrize('seq,pos,length,refr,alt,rwindow,awindow', [
    ('AACTAGCCTGCGGTCTGTGTTTCCCGACTTCTGAGTCATGGGGTTTCAATGCCTAT',
     5, 9, 'AGCCTGCGGT', 'A', 'ACTAGCCTGCGGTCTGT', 'ACTACTGT'),
    ('TTGAGATCGCGACGCTACTCTGAGCTCGGAGGAGCGGCATAAACGCGCCACCACCC',
     37, 4, 'GCATA', 'G', 'GCGGCATAAACG', 'GCGGAACG'),
    ('CCTTGGTGCCACGATCCGGCTATGGCGGAAGGGCACACCTAACCGCAACATTTGCC',
     14, 7, 'ATCCGGCT', 'A', 'ACGATCCGGCTATGG', 'ACGAATGG'),
    ('GGGTCCCAAGAGTCTGATTTCTAGCTTTTTATTTACACCCCGGTAGCAGGATCAGA',
     49, 5, 'GGATCA', 'G', 'GCAGGATCAGA', 'GCAGGA'),
])
def test_deletion(seq, pos, length, refr, alt, rwindow, awindow):
    testrefr, testalt, testrw, testaw = kevlar.gentrio.mutate_deletion(
        seq, pos, length, ksize=5
    )

    print('REFR', refr, testrefr, rwindow, testrw, file=sys.stderr)
    print('ALT', alt, testalt, awindow, testaw, file=sys.stderr)

    assert refr == testrefr
    assert alt == testalt
    assert rwindow == testrw
    assert awindow == testaw


def test_gen_muts():
    seqstream = kevlar.open(data_file('100kbx3.fa.gz'), 'r')
    sequences = kevlar.seqio.parse_seq_dict(seqstream)
    mutator = kevlar.gentrio.generate_mutations(sequences, seed=42)
    mutations = list(mutator)

    refrs = [m[2] for m in mutations]
    alts = [m[3] for m in mutations]

    print('DEBUG refrs', refrs, file=sys.stderr)
    print('DEBUG alts', alts, file=sys.stderr)

    testrefrs = [
        'T', 'A', 'GAGAGGGTTACATACGCAGAAAAAGACACACGCTACTGCCCCGCATAGCC', 'A',
        'C', 'A', 'A', 'C', 'CCAATTAACGTGGAAGCTCGGGGTACCAGCGGACTTAGTGTTCCCGTCT'
        'TCGGGTTTAATACTATACAGGCCGCAGATACGGCGAACGACGACCCCTTAATTAGACGTGAATTCACTT'
        'TAGCATCCTTGCTCCTTCCGAAAGACCCCTATGGAAGCAACTCCCGAAACCCCGAAAGGT', 'C'
    ]
    testalts = [
        'C', 'C', 'G', 'C', 'G', 'ACGCGGCATTACTTTCCCTTTTAATCGCATCCAGTAGTTACGGT'
        'GTAAGAGCTCGAGCGATTTTAGGTGTCAAAGGGAATTATTGAGGAAGTCTGGTATGCTGGGATAACCTG'
        'TACATCGACGGAATCTACCACCCTACACGGCCTAATACTCCCCGTTCCGTATACTGTATAGGACCTATA'
        'ATGGTTACTGACTGAACGAGTCCTACACTTCATGTCGCTGTAACATGGCGGATGTCCCAGTTGTTGTGC'
        'CGATATTCCCAATCTGTAAGTGTTCTATTCCGGC', 'C', 'T', 'C', 'T'
    ]

    assert refrs == testrefrs
    assert alts == testalts


def test_gen_with_inversions():
    with pytest.raises(NotImplementedError):
        list(kevlar.gentrio.generate_mutations({'1': 'ACGT'}, inversions=True))
