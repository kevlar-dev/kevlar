#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from tempfile import NamedTemporaryFile
import khmer
import kevlar
from kevlar.tests import data_file
from kevlar.simlike import get_abundances, abund_log_prob
from kevlar.simlike import likelihood_denovo
from kevlar.simlike import likelihood_false
from kevlar.simlike import likelihood_inherited
import pytest


@pytest.fixture
def minitrio():
    kid = khmer.Counttable(31, 1e6, 4)
    mom = khmer.Counttable(31, 1e6, 4)
    dad = khmer.Counttable(31, 1e6, 4)
    ref = khmer.SmallCounttable(31, 125000, 4)
    kid.consume_seqfile(data_file('minitrio/trio-proband.fq.gz'))
    mom.consume_seqfile(data_file('minitrio/trio-mother.fq.gz'))
    dad.consume_seqfile(data_file('minitrio/trio-father.fq.gz'))
    ref.consume_seqfile(data_file('minitrio/refr.fa'))
    return kid, mom, dad, ref


@pytest.fixture
def miniabund(minitrio):
    kid, mom, dad, ref = minitrio
    altseq = 'TGTCTCCCTCCCCTCCACCCCCAGAAATGGGTTTTTGATAGTCTTCCAAAGTTAGGGTAGT'
    refseq = 'TGTCTCCCTCCCCTCCACCCCCAGAAATGGCTTTTTGATAGTCTTCCAAAGTTAGGGTAGT'
    altabund, refrabund, ndropped = get_abundances(
        altseq, refseq, kid, (mom, dad), ref
    )
    assert ndropped == 3
    return altabund, refrabund


def test_get_abundances(minitrio):
    kid, mom, dad, ref = minitrio
    altseq = 'TGTCTCCCTCCCCTCCACCCCCAGAAATGGGTTTTTGATAGTCTTCCAAAGTTAGGGTAGT'
    refseq = 'TGTCTCCCTCCCCTCCACCCCCAGAAATGGCTTTTTGATAGTCTTCCAAAGTTAGGGTAGT'
    altabund, refrabund, ndropped = get_abundances(
        altseq, refseq, kid, (mom, dad), ref
    )
    assert ndropped == 3
    assert altabund == [
        [7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 9, 8, 8, 9, 9, 9, 7, 7, 8, 8, 8, 7,
         7, 7, 7, 7, 7],
        [1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0],
    ]
    assert refrabund == [2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                         1, 1, 1, 1, 2, 1, 1, 1, 1, 1]

    refseq = 'TGTCTCCCTCCCCTCCACCCCCAGAAATGGGAAATTTTTGATAGTCTTCCAAAGTTAGGGTAGT'
    altabund, refrabund, ndropped = get_abundances(
        altseq, refseq, kid, (mom, dad), ref
    )
    assert ndropped == 3
    assert refrabund == [None] * len(altabund[0])


def test_abund_log_prob():
    assert abund_log_prob(0, 3, refrabund=1) == pytest.approx(-10.51967)
    assert abund_log_prob(0, 4, refrabund=1) == pytest.approx(-14.02623)
    assert abund_log_prob(0, 4, refrabund=6) == pytest.approx(-6.85919)
    assert abund_log_prob(0, 4, refrabund=15) == pytest.approx(-3.19403)

    assert abund_log_prob(1, 1) == pytest.approx(-8.43023)
    assert abund_log_prob(1, 10) == pytest.approx(-3.08648)
    assert abund_log_prob(1, 15) == pytest.approx(-2.305233)
    assert abund_log_prob(1, 20) == pytest.approx(-3.08648)
    assert abund_log_prob(1, 10, mean=50.0, sd=9.9) == pytest.approx(-7.10969)
    assert abund_log_prob(1, 20, mean=50.0, sd=9.9) == pytest.approx(-3.02848)

    assert abund_log_prob(2, 1) == pytest.approx(-9.5687)
    assert abund_log_prob(2, 10) == pytest.approx(-6.12338)
    assert abund_log_prob(2, 30) == pytest.approx(-2.99838)
    assert abund_log_prob(2, 53) == pytest.approx(-7.13119)
    assert abund_log_prob(2, 29, mean=47.0, sd=9.3) == pytest.approx(-5.0220)
    assert abund_log_prob(2, 37, mean=47.0, sd=9.3) == pytest.approx(-3.727054)
    assert abund_log_prob(2, 43, mean=47.0, sd=9.3) == pytest.approx(-3.241449)


def test_likelihood_denovo(miniabund):
    altabund, refrabund = miniabund
    assert likelihood_denovo(altabund, refrabund) == pytest.approx(-221.90817)


def test_likelihood_false(miniabund):
    altabund, refrabund = miniabund
    assert likelihood_false(altabund, refrabund) == pytest.approx(-785.71390)


def test_likelihood_inherited(miniabund):
    altabund, refrabund = miniabund
    assert likelihood_inherited(altabund) == pytest.approx(-436.01119)


def test_joinlist():
    assert kevlar.simlike.joinlist([1, 2, 3, 4, 5]) == '1,2,3,4,5'
    assert kevlar.simlike.joinlist([]) == '.'


def test_simlike_main(minitrio):
    kid, mom, dad, ref = minitrio
    instream = kevlar.open(data_file('minitrio/calls.vcf'), 'r')
    reader = kevlar.vcf.VCFReader(instream)
    calculator = kevlar.simlike.simlike(
        reader, kid, (mom, dad), ref, samplelabels=('Kid', 'Mom', 'Dad')
    )
    calls = list(calculator)
    assert len(calls) == 1
    call = calls[0]
    assert float(call.attribute('LLDN')) == pytest.approx(-221.90817)
    assert call.format('Kid', 'ALTABUND') == ('7,6,6,6,6,6,6,6,6,6,7,9,8,8,9,9'
                                              ',9,7,7,8,8,8,7,7,7,7,7,7')


def test_simlike_main_no_labels(minitrio):
    kid, mom, dad, ref = minitrio
    instream = kevlar.open(data_file('minitrio/calls.vcf'), 'r')
    reader = kevlar.vcf.VCFReader(instream)
    calculator = kevlar.simlike.simlike(reader, kid, (mom, dad), ref)
    calls = list(calculator)
    assert len(calls) == 1
    labels = set(calls[0]._sample_data.keys())
    assert labels == set(('Case', 'Control1', 'Control2'))


@pytest.mark.parametrize('fmtstr,sampleargs', [
    (
        'FORMAT\tProband\tMother\tFather\n',
        ['--sample-labels', 'Proband', 'Mother', 'Father'],
    ),
    (
        'FORMAT\tCase\tControl1\tControl2\n',
        [],
    ),
])
def test_simlike_cli(fmtstr, sampleargs, minitrio, capsys):
    kid, mom, dad, ref = minitrio
    with NamedTemporaryFile(suffix='.ct') as kidct, \
            NamedTemporaryFile(suffix='.ct') as momct, \
            NamedTemporaryFile(suffix='.ct') as dadct, \
            NamedTemporaryFile(suffix='.sct') as refrsct:
        kid.save(kidct.name)
        mom.save(momct.name)
        dad.save(dadct.name)
        ref.save(refrsct.name)

        arglist = [
            'simlike', '--case', kidct.name,
            '--controls', momct.name, dadct.name, *sampleargs,
            '--refr', refrsct.name, data_file('minitrio/calls.vcf')
        ]
        print(arglist)
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.simlike.main(args)

    out, err = capsys.readouterr()
    assert fmtstr in out
    assert 'LIKESCORE=214.103' in out
    assert 'LLDN=-221.908;LLFP=-785.714;LLIH=-436.011' in out


def test_simlike_cli_bad_labels():
    arglist = [
        'simlike', '--case', 'kid.ct',
        '--controls', 'mom.ct', 'dad.ct',
        '--sample-labels', 'Proband', 'Mother', 'Father', 'Sibling',
        '--refr', 'refr.sct', data_file('minitrio/calls.vcf')
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    with pytest.raises(kevlar.simlike.KevlarSampleLabelingError) as sle:
        kevlar.simlike.main(args)
    assert 'provided 4 labels but 3 samples' in str(sle)
