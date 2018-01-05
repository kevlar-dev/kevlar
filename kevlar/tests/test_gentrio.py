#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from io import StringIO
import os
import pytest
import random
from tempfile import mkdtemp
from shutil import rmtree
import sys
import kevlar
from kevlar import MutableString
from kevlar.tests import data_file


def test_rng():
    rng = random.Random(1776)
    draws = [rng.randint(1, 100) for _ in range(5)]
    assert draws == [80, 75, 14, 46, 21]


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


def test_insertion_rng():
    seq = 'ATGCCTATAGATTCAGTAGTTACCAGAGGCAGTGGTGTTTGCCACGCCATTTCTACGCGA'
    rng = random.Random(2018)
    refr, alt, refrwindow, altwindow = kevlar.gentrio.mutate_insertion(
        seq, position=19, length=5, duplpos=44, rng=rng, ksize=11
    )
    assert refr == 'G'
    assert alt == 'GCCCCA'
    assert refrwindow == 'GATTCAGTAGTTACCAGAGG'
    assert altwindow == 'GATTCAGTAGCCCCATTACCAGAGG'


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
    mutator = kevlar.gentrio.generate_mutations(sequences, rng=42)
    mutations = list(mutator)

    refrs = [m._refr for m in mutations]
    alts = [m._alt for m in mutations]

    print('DEBUG refrs', refrs, file=sys.stderr)
    print('DEBUG alts', alts, file=sys.stderr)

    testrefrs = [
        'T', 'A', 'GAGAGGGTTACATACGCAGAAAAAGACACACGCTACTGCCCCGCATAGCC', 'A',
        'C', 'A', 'AGTTTGTCGACCTCGAACTTGCCGCC', 'A', 'T', 'C'
    ]
    testalts = [
        'C', 'C', 'G', 'C', 'G', 'ACGCTGCATTACTTTCCCTTTTAATGGCATCCAGTAGTTACGGT'
        'GTAAGAGCTCGAGCGATTTTAGGTGTCAAAGGGAATTATTGAGGAAGTATGGTATGCTGGGATGACCTG'
        'TACATCGACGGAATCTACCACCCTACTCGGCCTAATACTCCCCGTTCCGTATACTGTATAGGACCTATA'
        'ATGGTTACTGACTGAACGAGTCCTACACTTCATGTCGCTGTAACATGGCGGATGTCCCAGTTGTTGTGC'
        'AGATATTCCCAATCTGTAAGAGTTCTCTTCTGGC', 'A', 'T', 'TTTTCCTTGAGCTGTCTAGTG'
        'CCGTAGGTGAACCTCGTGGATATTAAGGCGTACGGACCATAGAGAGCTACCTCTAAGAGGTCCGGGGAT'
        'AGGTTCCGGTCGTTGGTACACCCTGTATACCCCGCGACAGTTTTGAGGTCACAGCGTAGAGGGTGTAAC'
        'GCAGTATGCAGGATTGGGATGGATGTGGACTTTGTATAGGTAGTATGTGTTTTTG', 'CAAATCGCCT'
        'TTACCTAATCGACAAATAGGTTCCGTTGTACTTGCTCTAACCGTCGGTGATTTGAACACTTTCCATGTT'
        'ACCACACTAG'
    ]

    assert refrs == testrefrs
    assert alts == testalts
    assert mutations[0].genotypes is None


def test_gen_with_inversions():
    with pytest.raises(NotImplementedError):
        list(kevlar.gentrio.generate_mutations({'1': 'ACGT'}, inversions=True))


def test_sim_var_geno_smoketest():
    seqstream = kevlar.open(data_file('100kbx3.fa.gz'), 'r')
    sequences = kevlar.seqio.parse_seq_dict(seqstream)
    ninh = random.randint(1, 10)
    ndenovo = random.randint(1, 10)
    simulator = kevlar.gentrio.simulate_variant_genotypes(
        sequences, ninh=ninh, ndenovo=ndenovo
    )
    variants = list(simulator)
    assert len(variants) == ninh + ndenovo


def test_sim_var_geno():
    seqstream = kevlar.open(data_file('100kbx3.fa.gz'), 'r')
    sequences = kevlar.seqio.parse_seq_dict(seqstream)
    simulator = kevlar.gentrio.simulate_variant_genotypes(
        sequences, ninh=2, ndenovo=2, rng=112358 ^ 853211
    )

    variants = list(simulator)
    seqids = [v.seqid for v in variants]
    positions = [v.position for v in variants]
    genotypes = [v.genotypes for v in variants]

    print('DEBUG', seqids, positions, genotypes)

    assert len(variants) == 4
    assert seqids == ['scaf3', 'scaf3', 'scaf3', 'scaf1']
    assert positions == [4936, 57391, 86540, 68352]
    assert genotypes == [
        ('0/1', '0/1', '1/0'),
        ('1/0', '0/1', '0/0'),
        ('0/1', '0/0', '0/0'),
        ('0/1', '0/0', '0/0')
    ]


def test_apply_mut_snv():
    contig = MutableString('ACGTACGTACGT')

    kevlar.gentrio.apply_mutation(contig, 5, 'C', 'G')
    assert contig == 'ACGTAGGTACGT'

    kevlar.gentrio.apply_mutation(contig, 5, 'G', 'A')
    assert contig == 'ACGTAAGTACGT'

    kevlar.gentrio.apply_mutation(contig, 0, 'A', 'T')
    assert contig == 'TCGTAAGTACGT'


def test_apply_mut_ins():
    contig = MutableString('ACGTACGTACGT')
    kevlar.gentrio.apply_mutation(contig, 5, 'A', 'AAAAA')
    assert contig == 'ACGTAAAAACGTACGT'

    contig = MutableString('CTTGAGACTTAGTAAAACCGTC')
    kevlar.gentrio.apply_mutation(contig, 7, 'A', 'ATTCTTGTT')
    assert contig == 'CTTGAGATTCTTGTTCTTAGTAAAACCGTC'


def test_apply_mut_del():
    contig = MutableString('ACGTACGTACGT')
    kevlar.gentrio.apply_mutation(contig, 5, 'ACGTAC', 'A')
    assert contig == 'ACGTAGT'


def test_gentrio_smoketest():
    seqstream = kevlar.open(data_file('100kbx3.fa.gz'), 'r')
    sequences = kevlar.seqio.parse_seq_dict(seqstream)
    outstreams = [StringIO(), StringIO(), StringIO()]
    mutator = kevlar.gentrio.gentrio(sequences, outstreams, ninh=2, ndenovo=1,
                                     seed=1985)
    variants = list(mutator)

    for variant in variants:
        print(variant.vcf, file=sys.stderr)

    for i in range(3):
        outstreams[i].seek(0)
    probandseqs = kevlar.seqio.parse_seq_dict(outstreams[0])
    motherseqs = kevlar.seqio.parse_seq_dict(outstreams[1])
    fatherseqs = kevlar.seqio.parse_seq_dict(outstreams[2])

    print(probandseqs['scaf1_haplo1'][variants[0].position])
    print(probandseqs['scaf1_haplo2'][variants[0].position])
    assert variants[0].genotypes[0] == '0/1'
    assert variants[0].refrwindow in probandseqs['scaf1_haplo1']
    assert variants[0].refrwindow not in probandseqs['scaf1_haplo2']
    assert variants[0].window not in probandseqs['scaf1_haplo1']
    assert variants[0].window in probandseqs['scaf1_haplo2']

    print(probandseqs['scaf3_haplo1'][variants[2].position])
    print(probandseqs['scaf3_haplo2'][variants[2].position])
    print(motherseqs['scaf3_haplo1'][variants[2].position])
    print(motherseqs['scaf3_haplo2'][variants[2].position])
    print(fatherseqs['scaf3_haplo1'][variants[2].position])
    print(fatherseqs['scaf3_haplo2'][variants[2].position])
    assert variants[2].window in probandseqs['scaf3_haplo1']
    assert variants[2].refrwindow in probandseqs['scaf3_haplo2']
    assert variants[2].refrwindow in motherseqs['scaf3_haplo1']
    assert variants[2].refrwindow in motherseqs['scaf3_haplo2']
    assert variants[2].refrwindow in fatherseqs['scaf3_haplo1']
    assert variants[2].window in fatherseqs['scaf3_haplo2']


def test_gentrio_cli_smoketest():
    tempdir = mkdtemp()
    prefix = os.path.join(tempdir, 'outfile')
    arglist = [
        'gentrio', '--prefix', prefix, data_file('100kbx3.fa.gz')
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.gentrio.main(args)
    rmtree(tempdir)


def test_gentrio_cli(capsys):
    tempdir = mkdtemp()
    prefix = os.path.join(tempdir, 'outfile')
    arglist = [
        'gentrio', '--prefix', prefix, '--vcf', '-', data_file('100kbx3.fa.gz')
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.gentrio.main(args)
    rmtree(tempdir)
    out, err = capsys.readouterr()

    outlines = out.strip().split('\n')
    numvariants = sum([1 for line in outlines if not line.startswith('#')])
    assert numvariants == 30


def test_mutagenize():
    rng = random.Random(1123581321)
    mutseq = kevlar.gentrio.mutagenize('GATTACA' * 3, rng, rate=0.1)
    assert mutseq == 'GATAACAGATTACAGATTACG'
