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
from kevlar.simlike import spanning_kmer_abundances, abund_log_prob
from kevlar.simlike import likelihood_denovo
from kevlar.simlike import likelihood_false
from kevlar.simlike import likelihood_inherited
import pytest


@pytest.fixture(scope="module")
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


@pytest.fixture(scope="module")
def miniabund(minitrio):
    kid, mom, dad, ref = minitrio
    altseq = 'TGTCTCCCTCCCCTCCACCCCCAGAAATGGGTTTTTGATAGTCTTCCAAAGTTAGGGTAGT'
    refseq = 'TGTCTCCCTCCCCTCCACCCCCAGAAATGGCTTTTTGATAGTCTTCCAAAGTTAGGGTAGT'
    altabund, refrabund, ndropped = spanning_kmer_abundances(
        altseq, refseq, kid, (mom, dad), ref
    )
    assert ndropped == 3
    return altabund, refrabund


@pytest.fixture(scope="module")
def caselowsketches():
    kid = kevlar.sketch.load(data_file('case-low-abund/kid.ct'))
    mom = kevlar.sketch.load(data_file('case-low-abund/mom.ct'))
    dad = kevlar.sketch.load(data_file('case-low-abund/dad.ct'))
    refr = kevlar.sketch.load(data_file('case-low-abund/refr.sct'))
    return kid, mom, dad, refr


@pytest.fixture(scope="module")
def ctrlhighsketches():
    kid = kevlar.sketch.load(data_file('ctrl-high-abund/cc57120.kid.sct'))
    mom = kevlar.sketch.load(data_file('ctrl-high-abund/cc57120.mom.sct'))
    dad = kevlar.sketch.load(data_file('ctrl-high-abund/cc57120.dad.sct'))
    refr = kevlar.sketch.load(data_file('ctrl-high-abund/cc57120.refr.sct'))
    return kid, mom, dad, refr


@pytest.fixture(scope="module")
def term_high_abund_trio():
    kid = kevlar.sketch.load(data_file('term-high-abund/proband.ct'))
    mom = kevlar.sketch.load(data_file('term-high-abund/mother.ct'))
    dad = kevlar.sketch.load(data_file('term-high-abund/father.ct'))
    refr = khmer.Nodetable(31, 1, 1)
    return kid, mom, dad, refr


def test_spanning_kmer_abundances(minitrio):
    kid, mom, dad, ref = minitrio
    altseq = 'TGTCTCCCTCCCCTCCACCCCCAGAAATGGGTTTTTGATAGTCTTCCAAAGTTAGGGTAGT'
    refseq = 'TGTCTCCCTCCCCTCCACCCCCAGAAATGGCTTTTTGATAGTCTTCCAAAGTTAGGGTAGT'
    altabund, refrabund, ndropped = spanning_kmer_abundances(
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
    altabund, refrabund, ndropped = spanning_kmer_abundances(
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
    assert likelihood_inherited(altabund) == pytest.approx(-438.31377)


def test_joinlist():
    assert kevlar.simlike.joinlist([1, 2, 3, 4, 5]) == '1,2,3,4,5'
    assert kevlar.simlike.joinlist([]) == '.'


def test_simlike_bad_windows(minitrio, capsys):
    kid, mom, dad, ref = minitrio
    instream = kevlar.open(data_file('minitrio/calls-badwindows.vcf'), 'r')
    reader = kevlar.vcf.VCFReader(instream)
    calculator = kevlar.simlike.simlike(
        reader, kid, (mom, dad), ref, samplelabels=('Kid', 'Mom', 'Dad')
    )
    calls = list(calculator)
    assert len(calls) == 5
    goodcalls = [c for c in calls if c.attribute('LIKESCORE') > float('-inf')]
    assert len(goodcalls) == 1
    assert len(goodcalls[0].window) == 61
    assert len(goodcalls[0].refrwindow) == 61

    out, err = capsys.readouterr()
    print('DEBUG', out)
    print('DEBUG', err)
    assert 'missing alt allele spanning window' in err
    assert 'missing refr allele spanning window' in err
    assert 'alt allele spanning window CCCCAGAAATGGGTTTTTGATAGTC' in err
    assert 'ref allele spanning window CCCCAGAAATGGCTTTTTGATAGTC' in err


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
    print(out)
    assert fmtstr in out
    assert 'LIKESCORE=216.406' in out
    assert 'LLDN=-221.908;LLFP=-785.714;LLIH=-438.314' in out


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


def test_simlike_fastmode():
    kid = kevlar.sketch.load(data_file('simlike-fast-mode/cc27.kid.ct'))
    mom = kevlar.sketch.load(data_file('simlike-fast-mode/cc27.mom.ct'))
    dad = kevlar.sketch.load(data_file('simlike-fast-mode/cc27.dad.ct'))
    refr = kevlar.sketch.load(data_file('simlike-fast-mode/cc27.refr.sct'))

    vcfin = kevlar.open(data_file('simlike-fast-mode/cc27.calls.vcf'), 'r')
    prelimcalls = kevlar.vcf.VCFReader(vcfin)
    scorer = kevlar.simlike.simlike(
        prelimcalls, kid, [mom, dad], refr, fastmode=True,
        samplelabels=['Proband', 'Mother', 'Father'],
    )
    calls = list(scorer)
    assert len(calls) == 4
    proband_abunds = [c.format('Proband', 'ALTABUND') for c in calls]
    filters = [c.filterstr for c in calls]
    assert proband_abunds == [None] * 4
    assert filters == [
        'LikelihoodFail;PassengerVariant', 'ControlAbundance;LikelihoodFail',
        'ControlAbundance;LikelihoodFail', 'LikelihoodFail;UserFilter'
    ]


@pytest.mark.parametrize('threshold,filterstatus', [
    (-10, False),
    (-1, False),
    (0, False),
    (None, False),
    (False, False),
    (1, True),
    (3, True),
    (5, True),
    (15, False),
    (49, False),
])
def test_simlike_ctrl_high_abund(threshold, filterstatus, ctrlhighsketches):
    kid, mom, dad, refr = ctrlhighsketches
    vcfin = kevlar.open(data_file('ctrl-high-abund/cc57120.calls.vcf'), 'r')
    prelimcalls = kevlar.vcf.VCFReader(vcfin)
    scorer = kevlar.simlike.simlike(
        prelimcalls, kid, [mom, dad], refr, samplelabels=['Kid', 'Mom', 'Dad'],
        ctrlabundhigh=threshold,
    )
    calls = list(scorer)
    assert len(calls) == 2
    print([c.filterstr for c in calls])
    for c in calls:
        assert ('ControlAbundance' in c.filterstr) is filterstatus


@pytest.mark.parametrize('casemin,abund,numfilt', [
    (6, -10, 0),
    (6, -1, 0),
    (6, 0, 0),
    (6, None, 0),
    (6, False, 0),
    (6, 5, 4),
    (7, 5, 5),
    (6, 4, 5),
    (6, 9, 4),
    (6, 10, 3),
])
def test_simlike_case_low_abund(casemin, abund, numfilt, caselowsketches):
    kid, mom, dad, refr = caselowsketches
    vcfin = kevlar.open(data_file('case-low-abund/calls.vcf.gz'), 'r')
    prelimcalls = kevlar.vcf.VCFReader(vcfin)
    scorer = kevlar.simlike.simlike(
        prelimcalls, kid, [mom, dad], refr, samplelabels=['Kid', 'Mom', 'Dad'],
        casemin=casemin, caseabundlow=abund,
    )
    calls = list(scorer)
    assert len(calls) == 5
    print([c.filterstr for c in calls])
    filtered = [c for c in calls if 'CaseAbundance' in c.filterstr]
    assert len(filtered) == numfilt


def test_simlike_min_like_score(ctrlhighsketches):
    kid, mom, dad, refr = ctrlhighsketches
    vcfin = kevlar.open(data_file('ctrl-high-abund/cc57120.calls.vcf'), 'r')
    prelimcalls = kevlar.vcf.VCFReader(vcfin)
    scorer = kevlar.simlike.simlike(
        prelimcalls, kid, [mom, dad], refr, samplelabels=['Kid', 'Mom', 'Dad'],
        ctrlabundhigh=0, caseabundlow=0, minlikescore=0.0,
    )
    calls = list(scorer)
    passing = [c for c in calls if c.filterstr == 'PASS']
    notpassing = [c for c in calls if c.filterstr != 'PASS']
    assert len(passing) == 1
    assert len(notpassing) == 1

    vcfin = kevlar.open(data_file('ctrl-high-abund/cc57120.calls.vcf'), 'r')
    prelimcalls = kevlar.vcf.VCFReader(vcfin)
    scorer = kevlar.simlike.simlike(
        prelimcalls, kid, [mom, dad], refr, samplelabels=['Kid', 'Mom', 'Dad'],
        ctrlabundhigh=0, caseabundlow=0, minlikescore=400.0,
    )
    calls = list(scorer)
    passing = [c for c in calls if c.filterstr == 'PASS']
    notpassing = [c for c in calls if c.filterstr != 'PASS']
    assert len(passing) == 0
    assert len(notpassing) == 2


@pytest.mark.parametrize('dodrop,filterstr', [
    (True, 'PASS'),
    (False, 'LikelihoodFail'),
])
def test_simlike_drop_outliers(dodrop, filterstr, term_high_abund_trio):
    kid, mom, dad, refr = term_high_abund_trio
    prelimcalls = kevlar.vcf.VCFReader(
        kevlar.open(data_file('term-high-abund/calls.vcf'), 'r')
    )
    scorer = kevlar.simlike.simlike(
        prelimcalls, kid, [mom, dad], refr, mu=30.0, sigma=10.0, casemin=5,
        ctrlmax=1, dropoutliers=dodrop, ambigthresh=0
    )
    for call in scorer:
        assert call.filterstr == filterstr


@pytest.mark.parametrize('ambigthresh,filterstr', [
    (64, 'PASS'),
    (0, 'PASS'),
    (10, 'AmbiguousCall'),
    (None, 'PASS'),
    (False, 'PASS'),
])
def test_simlike_ambig_threshold(ambigthresh, filterstr, term_high_abund_trio):
    kid, mom, dad, refr = term_high_abund_trio
    prelimcalls = kevlar.vcf.VCFReader(
        kevlar.open(data_file('term-high-abund/calls.vcf'), 'r')
    )
    scorer = kevlar.simlike.simlike(
        prelimcalls, kid, [mom, dad], refr, mu=30.0, sigma=10.0, casemin=5,
        ctrlmax=1, dropoutliers=True, ambigthresh=ambigthresh
    )
    calls = list(scorer)
    testcalls = [c for c in calls if c.attribute('PART') == '869']
    for call in testcalls:
        assert call.filterstr == filterstr
