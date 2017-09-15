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
import khmer
import kevlar
from kevlar.tests import data_file


def test_align():
    target = ('TAAATAAATATCTGGTGTTTGAGGCAAAAAGGCAGACTTAAATTCTAAATCACACCTGTGCTT'
              'CCAGCACTACCTTCAAGCGCAGGTTCGAGCCAGTCAGGCAGGGTACATAAGAGTCCATTGTGC'
              'CTGTATTATTTTGAGCAATGGCTAAAGTACCTTCACCCTTGCTCACTGCTCCCCCACTTCCTC'
              'AAGTCTCATCGTGTTTTTTTTAGAGCTAGTTTCTTAGTCTCATTAGGCTTCAGTCACCAT')
    query = ('TCTGGTGTTTGAGGCAAAAAGGCAGACTTAAATTCTAAATCACACCTGTGCTTCCAGCACTACC'
             'TTCAAGCGCAGGTTCGAGCCAGTCAGGACTGCTCCCCCACTTCCTCAAGTCTCATCGTGTTTTT'
             'TTTAGAGCTAGTTTCTTAGTCTCATTAGGCTTCAGTCACCATCATTTCTTATAGGAATACCA')
    assert kevlar.align(target, query) == '10D91M69D79M20I'


def test_call_62():
    qfile = data_file('ssc62.contig.augfasta')
    tfile = data_file('ssc62.gdna.fa')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    queryseqs = [record for record in qinstream]
    targetseqs = [record for record in khmer.ReadParser(tfile)]

    calls = [tup for tup in kevlar.call.call(targetseqs, queryseqs)]
    assert len(calls) == 1
    testcall = ('10_108283509-108283827', 'contig17;cc=1', '25D268M25D',
                '10:108283664:A->G')
    assert calls[0] == testcall


def test_call_106():
    qfile = data_file('ssc106.contig.augfasta.gz')
    tfile = data_file('ssc106.gdna.fa.gz')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    queryseqs = [record for record in qinstream]
    targetseqs = [record for record in khmer.ReadParser(tfile)]

    calls = [tup for tup in kevlar.call.call(targetseqs, queryseqs)]
    assert len(calls) == 1
    testcall = ('6_7464819-7465186', 'contig13;cc=1', '50D264M50D3M',
                '6:7464986:G->A')
    assert calls[0] == testcall


@pytest.mark.parametrize('targetfile,queryfile,cigar', [
    ('pico-7-refr.fa', 'pico-7-asmbl.fa', '10D83M190D75M20I1M'),
    ('pico-2-refr.fa', 'pico-2-asmbl.fa', '10D89M153I75M20I'),
])
def test_call_cli(targetfile, queryfile, cigar, capsys):
    target = data_file(targetfile)
    query = data_file(queryfile)
    args = kevlar.cli.parser().parse_args(['call', query, target])
    kevlar.call.main(args)

    out, err = capsys.readouterr()
    print(out.split('\n'))
    cigars = [line.split()[2] for line in out.strip().split('\n')]
    assert cigar in cigars
