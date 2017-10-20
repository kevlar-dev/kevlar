#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import re
import sys
import khmer
import kevlar
from kevlar.call import call, make_call
from kevlar.tests import data_file


def test_align():
    """Smoke test for ksw2 aligner"""
    target = ('TAAATAAATATCTGGTGTTTGAGGCAAAAAGGCAGACTTAAATTCTAAATCACACCTGTGCTT'
              'CCAGCACTACCTTCAAGCGCAGGTTCGAGCCAGTCAGGCAGGGTACATAAGAGTCCATTGTGC'
              'CTGTATTATTTTGAGCAATGGCTAAAGTACCTTCACCCTTGCTCACTGCTCCCCCACTTCCTC'
              'AAGTCTCATCGTGTTTTTTTTAGAGCTAGTTTCTTAGTCTCATTAGGCTTCAGTCACCAT')
    query = ('TCTGGTGTTTGAGGCAAAAAGGCAGACTTAAATTCTAAATCACACCTGTGCTTCCAGCACTACC'
             'TTCAAGCGCAGGTTCGAGCCAGTCAGGACTGCTCCCCCACTTCCTCAAGTCTCATCGTGTTTTT'
             'TTTAGAGCTAGTTTCTTAGTCTCATTAGGCTTCAGTCACCATCATTTCTTATAGGAATACCA')
    assert kevlar.align(target, query) == ('10D91M69D79M20I', 155)


@pytest.mark.parametrize('ccid,varcall', [
    ('5', 'seq1:185752:30D'),
    ('7', 'seq1:226611:190D'),
    ('9', 'seq1:1527139:I->TCCTGGTCTGCCACGGTTGACTTGCCTACATAT'),
])
def test_call_pico_indel(ccid, varcall):
    qfile = data_file('pico' + ccid + '.contig.augfasta')
    tfile = data_file('pico' + ccid + '.gdna.fa')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    queryseqs = [record for record in qinstream]
    targetseqs = [record for record in khmer.ReadParser(tfile)]

    calls = list(call(targetseqs, queryseqs))
    assert len(calls) == 1
    assert str(calls[0]) == varcall


@pytest.mark.parametrize('ccid,cigar,varcall', [
    ('62', '25D268M25D', '10:108283664:A->G'),
    ('106', '50D264M50D3M', '6:7464986:G->A'),
    ('223', '50D268M50D1M', '5:42345359:C->G'),
])
def test_call_ssc_isolated_snv(ccid, cigar, varcall):
    """
    Ensure isolated SNVs are called correctly.

    SNVs that are well separated from other variants have a distinct alignment
    signature as reflected in the CIGAR string reported by ksw2. They are
    either of the form "delete-match-delete" or "delete-match-delete-match",
    where the second match is very short (and spurious).
    """
    qfile = data_file('ssc' + ccid + '.contig.augfasta')
    tfile = data_file('ssc' + ccid + '.gdna.fa')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    queryseqs = [record for record in qinstream]
    targetseqs = [record for record in khmer.ReadParser(tfile)]

    calls = list(call(targetseqs, queryseqs))
    assert len(calls) == 1
    assert str(calls[0]) == varcall


def test_call_ssc_1bpdel():
    """Test 1bp deletion"""
    qfile = data_file('ssc218.contig.augfasta')
    tfile = data_file('ssc218.gdna.fa')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    query = [record for record in qinstream][0]
    target = [record for record in khmer.ReadParser(tfile)][0]
    variants = make_call(target, query, '50D132M1D125M50D', 31)

    assert isinstance(variants, list)
    assert len(variants) == 1
    assert str(variants[0]) == '6:23230160:1D'


def test_call_ssc_two_proximal_snvs():
    """
    Test two proximal SNVs

    Currently this serves as a negative control for calling isolated SNVs, but
    distinguishing which (if any) of a set of proximal SNVs is novel will be
    supported soon, and this test will need to be updated.
    """
    qfile = data_file('ssc107.contig.augfasta.gz')
    tfile = data_file('ssc107.gdna.fa.gz')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    query = [record for record in qinstream][0]
    target = [record for record in khmer.ReadParser(tfile)][0]

    variants = make_call(target, query, '25D263M25D', 31)
    assert len(variants) == 2


@pytest.mark.parametrize('targetfile,queryfile,cigar', [
    ('pico-7-refr.fa', 'pico-7-asmbl.fa', '10D83M190D75M20I1M'),
    ('pico-2-refr.fa', 'pico-2-asmbl.fa', '10D89M153I75M20I'),
])
def test_call_cli(targetfile, queryfile, cigar, capsys):
    """Smoke test for `kevlar call` cli"""
    target = data_file(targetfile)
    query = data_file(queryfile)
    args = kevlar.cli.parser().parse_args(['call', query, target])
    kevlar.call.main(args)

    out, err = capsys.readouterr()
    print(out)
    cigars = list()
    for line in out.strip().split('\n'):
        cigarmatch = re.search('CG=([^;\n]+)', line)
        if cigarmatch:
            cigar = cigarmatch.group(1)
            cigars.append(cigar)
    assert cigar in cigars


def test_snv_obj():
    snv = kevlar.call.VariantSNV('scaffold42', 10773, 'A', 'G')
    assert str(snv) == 'scaffold42:10773:A->G'
    vcfvalues = ['scaffold42', '10774', '.', 'A', 'G', '.', 'PASS', '.']
    assert snv.vcf == '\t'.join(vcfvalues)
    assert snv.cigar is None

    snv2 = kevlar.call.VariantSNV('chr5', 500, 'T', 'G', CG='10D200M10D')
    assert snv2.cigar == '10D200M10D'
    assert snv2.window is None


def test_indel_obj():
    """
    Test indel objects

    The coordinate used to construct the object is 0-based, but includes the
    nucleotide shared by the reference and alternate alleles. The str() output
    coordinate is increased by 1 to account for this nucleotide, while the VCF
    output is increased by 1 to transform to a 1-based system where the shared
    nucleotide is the point of reference.
    """
    indel1 = kevlar.call.VariantIndel('chr3', 8998622, 'GATTACA', 'G')
    assert str(indel1) == 'chr3:8998623:6D'
    vcfvalues = ['chr3', '8998623', '.', 'GATTACA', 'G', '.', 'PASS', '.']
    assert indel1.vcf == '\t'.join(vcfvalues)

    indel2 = kevlar.call.VariantIndel('chr6', 75522411, 'G', 'GATTACA')
    assert str(indel2) == 'chr6:75522412:I->ATTACA'
    vcfvalues = ['chr6', '75522412', '.', 'G', 'GATTACA', '.', 'PASS', '.']
    assert indel2.vcf == '\t'.join(vcfvalues)


def test_variant_kmers():
    #            variant here---------------|
    window = 'TTATTTTTAACAAAGGAGCAAAGGAGCAAAGGGCAAATACAATGAGGCAAAGATAGTCTCT'

    qfile = data_file('ssc223.contig.augfasta')
    tfile = data_file('ssc223.gdna.fa')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    queryseqs = [record for record in qinstream]
    targetseqs = [record for record in khmer.ReadParser(tfile)]

    calls = list(call(targetseqs, queryseqs))
    assert len(calls) == 1
    assert calls[0].window == window


@pytest.mark.parametrize('prefix,cigar,refrwindow,altwindow', [
    ('phony-snv-01', '25D98M25D',
        'GGGGGTGTCTGCGACCACAGCTGAACATGACGAAACGGGTG',
        'GGGGGTGTCTGCGACCACAGGTGAACATGACGAAACGGGTG'),
    ('phony-snv-02', '24D99M25D',
        'ATTCGTATTACCCCTGGGATTTGGGAGCTGGTCTATATAGG',
        'ATTCGTATTACCCCTGGGATATGGGAGCTGGTCTATATAGG'),
    ('phony-deletion-01', '25D28M8D49M25D',
        'GGCTCAAGACTAAAAAGACTGAGACTCGTTTTTGGTGACAAGCAGGGC',
        'GGCTCAAGACTAAAAAGACTTTTTTGGTGACAAGCAGGGC'),
    ('phony-deletion-02', '40D29M3D36M40D',
        'CATCATCTCGTAGGTTTGTCTAGTGCAAACAGAGTCCCCCTGC',
        'CATCATCTCGTAGGTTTGTCTGCAAACAGAGTCCCCCTGC'),
    ('phony-insertion-01', '10D34M7I49M10D1M',
        'CATCTGTTTTTCTCGAACTCGTATATTATCTATAAATTCC',
        'CATCTGTTTTTCTCGAACTCGATTACAGTATATTATCTATAAATTCC'),
    ('phony-insertion-02', '10D33M27I95M10D',
        'GCCAGGAAGTTTACGATAAGGTGTTGCCATTCGAAATGAC',
        'GCCAGGAAGTTTACGATAAGTATATATATATATATATATATATATATGTGTTGCCATTCGAAATGAC'),
])
def test_variant_window(prefix, cigar, refrwindow, altwindow):
    qfile = data_file(prefix + '.contig.fa')
    tfile = data_file(prefix + '.gdna.fa')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    query = [record for record in qinstream][0]
    target = [record for record in khmer.ReadParser(tfile)][0]

    variants = make_call(target, query, cigar, 21)
    assert len(variants) == 1
    assert variants[0].window == altwindow
    assert variants[0].refrwindow == refrwindow


def test_nocall():
    # Intentionally mismatched
    qfile = data_file('phony-deletion-01.contig.fa')
    tfile = data_file('phony-insertion-01.gdna.fa')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    query = [record for record in qinstream][0]
    target = [record for record in khmer.ReadParser(tfile)][0]

    variants = make_call(target, query, '25D5M22I5M46D8M13D2M35I', 21)
    assert len(variants) == 1
    assert variants[0].vcf == (
        'yourchr\t801\t.\t.\t.\t.\t.\t'
        'CG=25D5M22I5M46D8M13D2M35I;NC=inscrutablecigar;QN=contig4:cc=1;'
        'QS=AACTGGTGGGCTCAAGACTAAAAAGACTTTTTTGGTGACAAGCAGGGCGGCCTGCCCTTCCTGTAG'
        'TGCAAGAAAAT'
    )
