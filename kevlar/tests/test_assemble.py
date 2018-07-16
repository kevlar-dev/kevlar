#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import screed
from networkx import connected_components
import kevlar
import kevlar.__main__
from kevlar.tests import data_file


def test_fml_asm():
    fh = kevlar.open(data_file('reads2chain.fq.gz'), 'r')
    reads = [s for s in kevlar.parse_augmented_fastx(fh)]
    assert len(reads) == 16
    contigs = [c for c in kevlar.assembly.fml_asm(reads)]
    assert len(contigs) == 1
    assert contigs[0] == ('AAAACAAAAACAAACAAACAAAAAAAACTTCCTCCATTGGCACACAATGCA'
                          'ACTGCTTCCCTGTCTTGTACATGTGGAGATGTGATAAAGTAACTTCAGTGA'
                          'CAGTCAAATGTACTGTTACCTCAAAAAGTGCGATGCTTTCTTGCATAATTC'
                          'CTATCAATGTTCTATTTCACATATGTGATACATTATAAAATACATTTATCT'
                          'TTCACAGAATTCATTCTAGAGGGAAAATATTAACATGTTAGT')


@pytest.mark.parametrize('cc', [139, 27, 278, 327, 379])
def test_assembly_edgeless(cc):
    filename = 'edgeless/cc{:d}.afq.gz'.format(cc)
    fh = kevlar.open(data_file(filename), 'r')
    reads = [r for r in kevlar.parse_augmented_fastx(fh)]
    contigs = [c for c in kevlar.assembly.fml_asm(reads)]
    assert len(contigs) == 0


@pytest.mark.parametrize('cc,contig', [
    (110, 'CTTTAAGAGCTGTAACACTCACTGCGAAGGTCTGAGGCTTCATTCCTGAAGTCAGTGTAGACCATGA'
          'ACCCACGAGGAGGAACGAACAACTCTGGGTGCGCCACCTTTAAGAGCTGTAACACGGCTGGGCGCGG'
          'TGGCTCACGCCTGTAATCCTGGCACTTTGGGAGGCCGAGATGGGTGGATCACCAGGTCAGGAGATCA'
          'TAACCATCCTGGCTAACACGGTGAAACCCCATCTCTACT'),
    (206, 'AAATTATTTATGTGTCTAACTTTGTTACTAACATATGATAACTTTGAGGACAGAAGCAAGTCCCAGT'
          'CAACATTCTATATCCAACTGTTACCACAGAGCAAATAATAGGTGCGTAAACTGTTTGTTGATTGAGT'
          'ATAGTACTCAGGTGAGAATAAATGGAGAATGAAATAAAAGTGATATTGATCTGGGAGTATACTACAG'
          'TTCCCCTATCCAGGCAGAAAGTATATAATGCTTCTACAATAAGGATTGCAAAGCTACCAAAAAGGAA'
          'AAATGAAAACGTTGTAA'),
    (231, 'ACTACCCAAAGTATGTATTACATACTGTACATAAAATATCAAAGTACCCAAAATGTGTATTATATAC'
          'TCATCATAAAATATCAAACTACCCAAAGTATGTTTTACATACTGTACATAAAATATCAAAGTACCCA'
          'AAATGTGTATTACATACTGTACATAAAATATCAAACTACCCACAGTATGTATTACATACTTTACATA'
          'AAATATCACAGTACG'),
    (322, 'TCAGTATTTTGAACTGTAAAATGGGAAAAACAAAGCCAATACCACTTTTATCACTTATAAGTGATAT'
          'ATTTGTCTCTATTCATCTGTCTTCCTGCCTGTCTATAATAACATGGAGTATTTATATTTCTATGTTA'
          'GTAGTTAGCAATTAATAATTGCCCCATTAGATCTAATTAAATGAAGGAGCTTCTGCACAGCAAAAGA'
          'AACTATCATCGGAGTGAACAGGCAACCTACAGAATGGGAGAAAATTTTTGCAATCTACT'),
    (58, 'TAAAACAATAATTGCTAATATTCTTTAGGTAGCTGCTGTACAACAGCACTATGTTAAGAACTTCACAG'
         'GAATTGTCACATTCCCCATAAAACTTACATAATCCTACTATTATTTCCTGTTTCAGATAAGGAAAATG'
         'AAACCTCGCGAAGTTAGAAAACTTGTTCATTTTCATAGAGTTAATTAATCATTGGAACCAGGATATAA'
         'AGTCAAGGTGTGGGAATCTATACCTAGAGTGTAACCATGACATGCATCTCAAAAACCAACAATGGAAT'
         'CATAGAAGA'),
])
def test_assembly(cc, contig):
    filename = 'fml/cc{:d}.afq.gz'.format(cc)
    fh = kevlar.open(data_file(filename), 'r')
    reads = [r for r in kevlar.parse_augmented_fastx(fh)]
    contigs = [c for c in kevlar.assembly.fml_asm(reads)]
    assert len(contigs) == 1
    assert contigs[0] == contig


def test_assemble_main(capsys):
    cliargs = ['assemble', data_file('var1.reads.augfastq')]
    args = kevlar.cli.parser().parse_args(cliargs)
    kevlar.assemble.main(args)
    out, err = capsys.readouterr()
    contig = ('GTCCTTGAGTCCATTAGAGACGGCTTCCGCCGTAGGCCCACTTCCTTAAAGTCGAGACTTCTA'
              'AAAACCGGGGTGTAACTCTTTTATTACAAAGCGACTATCCACCTGTAAGGACAGTGATA')
    print('DEBUG', contig)
    print('DEBUG', out)
    assert contig in out or kevlar.revcom(contig) in out


def test_assemble_no_edges(capsys):
    cliargs = ['assemble', data_file('asmbl-no-edges.augfastq.gz')]
    args = kevlar.cli.parser().parse_args(cliargs)
    kevlar.assemble.main(args)
    out, err = capsys.readouterr()
    assert out == ''
