#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from sys import stderr
import pytest
import kevlar
from kevlar.mutate import Mutation
from kevlar.tests import data_file


def test_load_mutations_x():
    instream = kevlar.open(data_file('muts-x.txt'), 'r')
    mutations = kevlar.mutate.load_mutations(instream, stderr)
    assert len(mutations) == 1
    assert '1' in mutations
    assert len(mutations['1']) == 1

    mut = mutations['1'][0]
    assert mut == Mutation(seq='1', pos=441274, type='snv', data='3')


def test_load_mutations_y():
    instream = kevlar.open(data_file('muts-y.tsv'), 'r')
    mutations = kevlar.mutate.load_mutations(instream, stderr)
    assert len(mutations) == 3

    assert 'scaffold399' in mutations
    assert len(mutations['scaffold399']) == 1
    mut = mutations['scaffold399'][0]
    assert mut == Mutation(seq='scaffold399', pos=685357, type='ins',
                           data='AGCTACCCCAGTGAGTCGGTAATGTGATC')

    assert 'scaffold982' in mutations
    assert len(mutations['scaffold982']) == 1
    mut = mutations['scaffold982'][0]
    assert mut == Mutation(seq='scaffold982', pos=108754, type='del',
                           data='23')

    assert 'scaffold1102' in mutations
    assert len(mutations['scaffold1102']) == 1
    mut = mutations['scaffold1102'][0]
    assert mut == Mutation(seq='scaffold1102', pos=260686, type='snv',
                           data='1')


def test_load_mutations_z():
    instream = kevlar.open(data_file('muts-z.csv'), 'r')
    with pytest.raises(ValueError) as ve:
        mutations = kevlar.mutate.load_mutations(instream, stderr)
    assert 'error parsing mutation' in str(ve)


def test_mutate_snv():
    mutation = Mutation(seq='contig', pos=5, type='snv', data='1')
    contig = 'ACGTACGTACGT'
    assert kevlar.mutate.mutate_snv(contig, mutation) == 'ACGTAGGTACGT'

    mutation = Mutation(seq='contig', pos=5, type='snv', data='-1')
    assert kevlar.mutate.mutate_snv(contig, mutation) == 'ACGTAAGTACGT'

    mutation = Mutation(seq='contig', pos=0, type='snv', data='-1')
    assert kevlar.mutate.mutate_snv(contig, mutation) == 'TCGTACGTACGT'


def test_mutate_ins():
    mutation = Mutation(seq='contig', pos=5, type='ins', data='AAAA')
    contig = 'ACGTACGTACGT'
    mutcontig = 'ACGTAAAAACGTACGT'
    assert kevlar.mutate.mutate_insertion(contig, mutation) == mutcontig


def test_mutate_del():
    mutation = Mutation(seq='contig', pos=5, type='ins', data='5')
    contig = 'ACGTACGTACGT'
    assert kevlar.mutate.mutate_deletion(contig, mutation) == 'ACGTAGT'


def test_mutate_inv():
    mutation = Mutation(seq='contig', pos=5, type='inv', data='5')
    contig = 'ACGTACGTACGT'
    assert kevlar.mutate.mutate_inversion(contig, mutation) == 'ACGTACATGCGT'


def test_mutate_bogus():
    instream = kevlar.open(data_file('muts-w.txt'), 'r')
    with pytest.raises(ValueError) as ve:
        mutations = kevlar.mutate.load_mutations(instream, stderr)
    assert 'invalid variant type "slippage"' in str(ve)


def test_mutate_main(capsys):
    genome = data_file('mut-genome.fa')
    muts = data_file('mut-genome.txt')
    args = kevlar.cli.parser().parse_args(['mutate', muts, genome])
    kevlar.mutate.main(args)
    out, err = capsys.readouterr()

    contig1 = '>contig1\nGTACGGCTATTGTCTGAGCTCTTTTTAAGACTAATACGCGCTGGCTCACGGAA'
    assert contig1 in out

    contig2 = '>contig2\nGTCATGAACTGACTCGCACGCGCTTCGGAAATTGCCGTATGATATGAC'
    assert contig2 in out

    contig3 = '>contig3\nAGTCGAGTATTGTGGCATAAGCGGAACA'
    assert contig3 in out
