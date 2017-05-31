#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import kevlar
from kevlar import VariantSet


@pytest.fixture
def basicvset():
    vs = VariantSet()
    vs.add_kmer('TTTTT', 'read1')
    vs.add_kmer('TTTTT', 'read2')
    vs.add_kmer('TTTTT', 'read3')
    vs.add_kmer('TTTTA', 'read4')
    vs.add_contig('TTAAAAATT', 'TTTTT')
    vs.add_contig('TTAAAAATT', 'TTTTA')
    return vs


def test_vset_kmers(basicvset):
    assert basicvset.nkmers == 2
    assert basicvset.nreads == 4
    assert sorted(list(basicvset.kmers.keys())) == ['AAAAA', 'TAAAA']
    assert basicvset.kmers['AAAAA'] == set(['read1', 'read2', 'read3'])
    assert basicvset.kmers['TAAAA'] == set(['read4'])


def test_vset_contigs(basicvset):
    assert basicvset.ncontigs == 1
    assert list(basicvset.contigs.keys()) == ['AATTTTTAA']


def test_vset_write(basicvset, capsys):
    from sys import stdout
    basicvset.write(stdout)
    out, err = capsys.readouterr()
    outtest = ('Contig,ContigRevCom\tNumReads\tNumKmers\tReads\tKmers\n'
               'AATTTTTAA,TTAAAAATT\t4\t2\tread1,read2,read3,read4'
               '\tAAAAA,TAAAA\n')
    assert out == outtest
