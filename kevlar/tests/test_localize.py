#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from io import StringIO
import pytest
import kevlar
from kevlar.localize import Localizer, KevlarRefrSeqNotFoundError
from kevlar.localize import decompose_seeds, contigs_2_seeds, get_seed_matches
from kevlar.localize import cutout, localize
from kevlar.reference import ReferenceCutout
from kevlar.seqio import Record
from kevlar.tests import data_file
from tempfile import NamedTemporaryFile


def test_localizer_simple():
    intervals = Localizer(seedsize=25)
    assert list(intervals.get_cutouts()) == []

    intervals.add_seed_match('chr1', 100)
    intervals.add_seed_match('chr1', 115)
    intervals.add_seed_match('chr2', 200)
    intervals.add_seed_match('chr2', 205)
    intervals.add_seed_match('chr2', 207)
    intervals.add_seed_match('chr2', 235008)
    intervals.add_seed_match('chr2', 235075)
    testint = [c.interval for c in intervals.get_cutouts()]
    print('DEBUG', testint)
    assert testint == [
        ('chr1', 100, 140),
        ('chr2', 200, 232),
        ('chr2', 235008, 235100)
    ]


def test_localizer_incl_excl():
    intervals = Localizer(seedsize=25)
    intervals.add_seed_match('1', 100)
    intervals.add_seed_match('1', 120)
    intervals.add_seed_match('12', 200)
    intervals.add_seed_match('12', 209)
    intervals.add_seed_match('12', 213)
    intervals.add_seed_match('X', 1234)
    intervals.add_seed_match('X', 1245)
    intervals.add_seed_match('Un', 13579)
    intervals.add_seed_match('Un', 13597)

    testint = [c.interval for c in intervals.get_cutouts()]
    assert sorted(testint) == [
        ('1', 100, 145),
        ('12', 200, 238),
        ('Un', 13579, 13622),
        ('X', 1234, 1270),
    ]

    intervals.exclpattern = 'Un'
    testint = [c.interval for c in intervals.get_cutouts()]
    assert sorted(testint) == [
        ('1', 100, 145),
        ('12', 200, 238),
        ('X', 1234, 1270),
    ]

    intervals.inclpattern = r'^\d+$'
    testint = [c.interval for c in intervals.get_cutouts()]
    assert sorted(testint) == [
        ('1', 100, 145),
        ('12', 200, 238),
    ]


def test_get_cutouts_basic():
    intervals = Localizer(seedsize=10)
    intervals.add_seed_match('bogus-genome-chr2', 10)
    refrstream = open(data_file('bogus-genome/refr.fa'), 'r')
    seqs = kevlar.seqio.parse_seq_dict(refrstream)
    cutouts = list(intervals.get_cutouts(refrseqs=seqs))
    assert len(cutouts) == 1
    assert cutouts[0].defline == 'bogus-genome-chr2_10-20'
    assert cutouts[0].sequence == 'GTTACATTAC'


def test_get_cutouts_basic_2():
    intervals = Localizer(seedsize=21)
    intervals.add_seed_match('simple', 49)
    intervals.add_seed_match('simple', 52)
    intervals.add_seed_match('simple', 59)
    refrstream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    seqs = kevlar.seqio.parse_seq_dict(refrstream)
    cutouts = list(intervals.get_cutouts(refrseqs=seqs, delta=5))
    assert len(cutouts) == 1
    assert cutouts[0].defline == 'simple_44-85'
    assert cutouts[0].sequence == 'AATACTATGCCGATTTATTCTTACACAATTAAATTGCTAGT'


def test_get_cutouts_basic_3():
    intervals = Localizer(seedsize=21)
    intervals.add_seed_match('simple', 40)
    intervals.add_seed_match('simple', 80)
    intervals.add_seed_match('simple', 120)
    intervals.add_seed_match('simple', 500)
    refrstream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    seqs = kevlar.seqio.parse_seq_dict(refrstream)
    cutouts = list(
        intervals.get_cutouts(refrseqs=seqs, clusterdist=None, delta=10)
    )
    assert len(cutouts) == 1
    assert cutouts[0].defline == 'simple_30-531'
    assert len(cutouts[0].sequence) == 501


def test_get_cutouts_large_span():
    refrstream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    seqs = kevlar.seqio.parse_seq_dict(refrstream)

    intervals = Localizer(seedsize=21)
    intervals.add_seed_match('simple', 100)
    intervals.add_seed_match('simple', 200)

    cutouts = intervals.get_cutouts(refrseqs=seqs, clusterdist=50, delta=25)
    deflines = [c.defline for c in cutouts]
    assert deflines == ['simple_75-146', 'simple_175-246']

    cutouts = intervals.get_cutouts(refrseqs=seqs, clusterdist=100, delta=50)
    deflines = [c.defline for c in cutouts]
    assert deflines == ['simple_50-271']


def test_get_cutouts_missing_seq():
    intervals = Localizer(seedsize=21)
    intervals.add_seed_match('simple', 100)
    intervals.add_seed_match('simple', 200)
    intervals.add_seed_match('TheCakeIsALie', 42)
    intervals.add_seed_match('TheCakeIsALie', 100)
    intervals.add_seed_match('TheCakeIsALie', 77)
    refrstream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    seqs = kevlar.seqio.parse_seq_dict(refrstream)
    with pytest.raises(KevlarRefrSeqNotFoundError) as rnf:
        list(intervals.get_cutouts(refrseqs=seqs))
    assert 'TheCakeIsALie' in str(rnf)


def test_extract_regions_boundaries():
    refrstream = open(data_file('simple-genome-ctrl1.fa'), 'r')
    seqs = kevlar.seqio.parse_seq_dict(refrstream)

    intervals = Localizer(seedsize=31)
    intervals.add_seed_match('simple', 15)
    cutouts = list(intervals.get_cutouts(refrseqs=seqs, delta=20))
    assert len(cutouts) == 1
    assert cutouts[0].defline == 'simple_0-66'

    intervals = Localizer(seedsize=31)
    intervals.add_seed_match('simple', 925)
    intervals.add_seed_match('simple', 955)
    intervals.add_seed_match('simple', 978)
    cutouts = list(intervals.get_cutouts(refrseqs=seqs, delta=20))
    assert len(cutouts) == 1
    assert cutouts[0].defline == 'simple_905-1000'


@pytest.mark.parametrize('X,numtargets', [
    (100000, 1),
    (10000, 5),
    (1000, 33),
    (0, 1),
    (None, 33),
])
def test_maxdiff(X, numtargets):
    contigstream = kevlar.parse_partitioned_reads(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file('maxdiff-contig.augfasta'), 'r')
        )
    )
    refrfile = data_file('maxdiff-refr.fa.gz')
    targeter = kevlar.localize.localize(contigstream, refrfile, seedsize=51,
                                        delta=50, maxdiff=X)
    targets = [cutout for partid, cutout in targeter]
    print([t.defline for t in targets])
    assert len(targets) == numtargets


@pytest.mark.parametrize('incl,excl,output', [
    (None, None, '>seq1_10-191'),
    (r'seq1', None, '>seq1_10-191'),
    (None, 'seq1', 'WARNING: no reference matches'),
    (r'chr[XY]', None, 'WARNING: no reference matches'),
    (None, r'b0Gu$', '>seq1_10-191'),
])
def test_main(incl, excl, output, capsys):
    contig = data_file('localize-contig.fa')
    refr = data_file('localize-refr.fa')
    arglist = ['localize', '--seed-size', '23', '--delta', '50', contig, refr]
    args = kevlar.cli.parser().parse_args(arglist)
    args.include = incl
    args.exclude = excl
    kevlar.localize.main(args)
    out, err = capsys.readouterr()
    print(out)
    assert output in out or output in err


def test_main_no_matches(capsys):
    contig = data_file('localize-contig-bad.fa')
    refr = data_file('localize-refr.fa')
    arglist = ['localize', '--seed-size', '23', contig, refr]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.localize.main(args)
    out, err = capsys.readouterr()
    assert 'WARNING: no reference matches' in err


def test_decompose_seeds():
    assert list(decompose_seeds('GATTACA', 5)) == ['GATTA', 'ATTAC', 'TTACA']
    assert list(decompose_seeds('GATTACA', 3)) == ['GAT', 'ATT', 'TTA', 'TAC',
                                                   'ACA']


def test_contigs_2_seeds():
    seedfile = StringIO()
    partitions = [[Record(name='seq', sequence='GATTACA')]]
    contigs_2_seeds(partitions, seedfile, seedsize=5)
    testoutput = '>seed0\nATTAC\n>seed1\nGATTA\n>seed2\nTGTAA\n'
    assert seedfile.getvalue() == testoutput


def test_get_seed_matches():
    seedfasta = (
        '>seed0\nATCTGTTCTTGGCCAATAGAAAAAGCAAGGAGCCCTGAAAGACTCACAGTG\n'
        '>seed1\nAAAAGGAAATGTTAACAACAAAATCACACAGATAAACCATCACAAGATCTG\n'
        '>seed2\nGATTCTAGGAGCTTGTTACTGCTGCTGAAAAAGGAAATGTTAACAACAAAA\n'
        '>seed3\nAACCAATAGAGGTCCACAGAAGTATATATAATCTGTTCTTGGCCAATAGAA\n'
        '>seed4\nTTGTGTGTAAAAACCAATAGAGGTCCACAGAAGTATATATAATCTGTTCTT\n'
        '>seed5\nAAGATACTATAATATGTTTCCCTGAGCACACCCCTTCGAAAGAGCAGAATT\n'
    )
    with NamedTemporaryFile(suffix='.fa', mode='w') as seedfile:
        print(seedfasta, file=seedfile, flush=True)
        refrfile = data_file('fiveparts-refr.fa.gz')
        seed_matches = get_seed_matches(seedfile.name, refrfile, seedsize=51)
        print(seed_matches)
        assert seed_matches == {
            'AACCAATAGAGGTCCACAGAAGTATATATAATCTGTTCTTGGCCAATAGAA': {('seq1',
                                                                     284819)},
            'AAGATACTATAATATGTTTCCCTGAGCACACCCCTTCGAAAGAGCAGAATT': {('seq1',
                                                                     284722)},
            'ATCTGTTCTTGGCCAATAGAAAAAGCAAGGAGCCCTGAAAGACTCACAGTG': {('seq1',
                                                                     284849)},
            'AAGAACAGATTATATATACTTCTGTGGACCTCTATTGGTTTTTACACACAA': {('seq1',
                                                                     284808)},
        }


def test_get_seed_matches_no_matches():
    seedfasta = (
        '>seed0\nAAAAGGAAATGTTAACAACAAAATCACACAGATAAACCATCACAAGATCTG\n'
        '>seed1\nGATTCTAGGAGCTTGTTACTGCTGCTGAAAAAGGAAATGTTAACAACAAAA\n'
    )
    with NamedTemporaryFile(suffix='.fa', mode='w') as seedfile:
        print(seedfasta, file=seedfile, flush=True)
        refrfile = data_file('fiveparts-refr.fa.gz')
        seed_matches = get_seed_matches(seedfile.name, refrfile, seedsize=51)
        assert seed_matches == {}


def test_localize_new():
    refr_file = data_file('fiveparts-refr.fa.gz')
    contig_file = data_file('fiveparts.contigs.augfasta.gz')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contig_file, 'r'))
    pstream = kevlar.parse_partitioned_reads(contigstream)
    localizer = localize(pstream, refr_file, seedsize=51, debug=True)
    cutoutdata = list(localizer)
    partids = [partid for partid, gdna in cutoutdata]
    gdnas = [gdna for partid, gdna in cutoutdata]
    deflines = [g.defline for g in gdnas]
    assert partids == ['1', '1', '2', '3', '4', '5']
    assert sorted(deflines) == sorted([
        'seq1_284663-284950', 'seq1_1924681-1925049', 'seq1_1660589-1660884',
        'seq1_2315741-2316037', 'seq1_2321099-2321322', 'seq1_593102-593389'
    ])


def test_localize_no_match():
    log = StringIO()
    refr_file = data_file('fiveparts-refr.fa.gz')
    contig_file = data_file('wasp-pass.contig.augfasta')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contig_file, 'r'))
    pstream = kevlar.parse_partitioned_reads(contigstream)
    localizer = localize(pstream, refr_file, seedsize=41, debug=True,
                         logstream=log)
    cutoutdata = list(localizer)
    assert cutoutdata == []
    assert 'WARNING: no reference matches' in log.getvalue()


@pytest.mark.parametrize('partid,testdeflines', [
    ('1', ['seq1_1924681-1925049', 'seq1_284663-284950']),
    ('4', ['seq1_2321099-2321322'])
])
def test_localize_new_single_partition(partid, testdeflines):
    refr_file = data_file('fiveparts-refr.fa.gz')
    contig_file = data_file('fiveparts.contigs.augfasta.gz')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contig_file, 'r'))
    pstream = kevlar.parse_single_partition(contigstream, partid)
    localizer = localize(pstream, refr_file, seedsize=51)
    cutoutdata = list(localizer)
    partids = [partid for partid, gdna in cutoutdata]
    gdnas = [gdna for partid, gdna in cutoutdata]
    deflines = sorted([g.defline for g in gdnas])
    assert deflines == testdeflines


def test_localize_cli(capsys):
    refr_file = data_file('fiveparts-refr.fa.gz')
    contig_file = data_file('fiveparts.contigs.augfasta.gz')

    arglist = ['localize', '--part-id', '2', contig_file, refr_file]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.localize.main(args)
    out, err = capsys.readouterr()
    assert out == (
        '>seq1_1660589-1660884 kvcc=2\n'
        'GATAGATCTCCAAGAATTTTATACAGCAGGGCCCTGAGAATGAGCATGGAAGTGAATTTATTAGCCAGT'
        'GACAGTCACTTCACACTCTTCCTATATCAAAATTGAAGCCCAGGCTGGAGGTGGGCAGGGGTAGTACTT'
        'TTATGGACTGGACAGGGCGTAATCCCACCTGGGCGTGGGAGGAATATAAAAATAACCTTTAATTAATTC'
        'TGTCTGTAATTTATCTATGGGATGGGGTTGTTCAGAGAAGACTTCAATACCAGTTATTTAAGCCTGACC'
        'CTGGCTTGCCTTGACCCCA\n'
    )

    arglist = ['localize', contig_file, refr_file]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.localize.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    assert len(outlines) == 12
