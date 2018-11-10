#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from io import StringIO
import pytest
import kevlar
from kevlar.cutout import decompose_seeds, contigs_2_seeds, get_seed_matches
from kevlar.cutout import cutout
from kevlar.seqio import Record
from kevlar.tests import data_file
from tempfile import NamedTemporaryFile


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
        assert seed_matches == {
            'AACCAATAGAGGTCCACAGAAGTATATATAATCTGTTCTTGGCCAATAGAA': ('seq1',
                                                                    284819),
            'AAGATACTATAATATGTTTCCCTGAGCACACCCCTTCGAAAGAGCAGAATT': ('seq1',
                                                                    284722),
            'ATCTGTTCTTGGCCAATAGAAAAAGCAAGGAGCCCTGAAAGACTCACAGTG': ('seq1',
                                                                    284849),
            'TTGTGTGTAAAAACCAATAGAGGTCCACAGAAGTATATATAATCTGTTCTT': ('seq1',
                                                                    284808),
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


def test_cutout():
    refr_file = data_file('fiveparts-refr.fa.gz')
    contig_file = data_file('fiveparts.contigs.augfasta.gz')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contig_file, 'r'))
    pstream = kevlar.parse_partitioned_reads(contigstream)
    localizer = cutout(pstream, refr_file, seedsize=51, debug=True)
    cutoutdata = list(localizer)
    partids = [partid for partid, gdna in cutoutdata]
    gdnas = [gdna for partid, gdna in cutoutdata]
    deflines = [g.defline for g in gdnas]
    assert partids == ['1', '1', '2', '3', '4', '5']
    assert deflines == [
        'seq1_284664-284950', 'seq1_1924681-1925048', 'seq1_1660589-1660884',
        'seq1_2315741-2316037', 'seq1_2321101-2321322', 'seq1_593102-593386'
    ]


def test_cutout():
    log = StringIO()
    refr_file = data_file('fiveparts-refr.fa.gz')
    contig_file = data_file('wasp-pass.contig.augfasta')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contig_file, 'r'))
    pstream = kevlar.parse_partitioned_reads(contigstream)
    localizer = cutout(pstream, refr_file, seedsize=41, debug=True,
                       logstream=log)
    cutoutdata = list(localizer)
    assert cutoutdata == []
    assert 'WARNING: no reference matches' in log.getvalue()


@pytest.mark.parametrize('partid,testdeflines', [
    ('1', ['seq1_284664-284950', 'seq1_1924681-1925048']),
    ('4', ['seq1_2321101-2321322'])
])
def test_cutout_single_partition(partid, testdeflines):
    refr_file = data_file('fiveparts-refr.fa.gz')
    contig_file = data_file('fiveparts.contigs.augfasta.gz')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contig_file, 'r'))
    pstream = kevlar.parse_single_partition(contigstream, partid)
    localizer = cutout(pstream, refr_file, seedsize=51)
    cutoutdata = list(localizer)
    partids = [partid for partid, gdna in cutoutdata]
    gdnas = [gdna for partid, gdna in cutoutdata]
    deflines = [g.defline for g in gdnas]
    assert deflines == testdeflines


def test_cutout_cli(capsys):
    refr_file = data_file('fiveparts-refr.fa.gz')
    contig_file = data_file('fiveparts.contigs.augfasta.gz')

    arglist = ['cutout', '--part-id', '2', contig_file, refr_file]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.cutout.main(args)
    out, err = capsys.readouterr()
    assert out == (
        '>seq1_1660589-1660884 kvcc=2\n'
        'GATAGATCTCCAAGAATTTTATACAGCAGGGCCCTGAGAATGAGCATGGAAGTGAATTTATTAGCCAGT'
        'GACAGTCACTTCACACTCTTCCTATATCAAAATTGAAGCCCAGGCTGGAGGTGGGCAGGGGTAGTACTT'
        'TTATGGACTGGACAGGGCGTAATCCCACCTGGGCGTGGGAGGAATATAAAAATAACCTTTAATTAATTC'
        'TGTCTGTAATTTATCTATGGGATGGGGTTGTTCAGAGAAGACTTCAATACCAGTTATTTAAGCCTGACC'
        'CTGGCTTGCCTTGACCCCA\n'
    )

    arglist = ['cutout', contig_file, refr_file]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.cutout.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    assert len(outlines) == 12
