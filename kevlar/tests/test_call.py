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
import re
import sys
import khmer
import kevlar
from kevlar.call import call, alignment_interpretable, VariantMapping
from kevlar.tests import data_file
import screed


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


def test_cigar_scrutability():
    assert alignment_interpretable('6M73I13M4I5M4D63M11I') is False
    assert alignment_interpretable('29I17M7I3M24I36M10D7M11I14M23I4M') is False
    assert alignment_interpretable('57I16M7I3M8D33M1I28M27I3M') is False

    assert alignment_interpretable('25D100M25D') is True
    assert alignment_interpretable('50I50M50I50M24D1M') is True
    assert alignment_interpretable('47I127M30D1M') is True
    assert alignment_interpretable('30D100M16D67M30D') is True


@pytest.mark.parametrize('ccid,varcall', [
    ('5', 'seq1:185752:30D'),
    ('7', 'seq1:226611:190D'),
    ('9', 'seq1:1527139:I->TCCTGGTCTGCCACGGTTGACTTGCCTACATAT'),
])
def test_call_pico_indel(ccid, varcall):
    qfile = data_file('pico' + ccid + '.contig.augfasta')
    tfile = data_file('pico' + ccid + '.gdna.fa')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    queries = [record for record in qinstream]
    tinstream = kevlar.reference.load_refr_cutouts(kevlar.open(tfile, 'r'))
    targets = [record for record in tinstream]

    calls = list(call(targets, queries))
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
    tinstream = kevlar.reference.load_refr_cutouts(kevlar.open(tfile, 'r'))
    targetseqs = [record for record in tinstream]

    calls = list(call(targetseqs, queryseqs))
    assert len(calls) == 1
    assert str(calls[0]) == varcall


def test_call_ssc_1bpdel():
    """Test 1bp deletion"""
    qfile = data_file('ssc218.contig.augfasta')
    tfile = data_file('ssc218.gdna.fa')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    query = [record for record in qinstream][0]
    tinstream = kevlar.reference.load_refr_cutouts(kevlar.open(tfile, 'r'))
    target = [record for record in tinstream][0]
    aln = VariantMapping(query, target, 1e6, '50D132M1D125M50D')
    variants = aln.call_variants(31)

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
    tinstream = kevlar.reference.load_refr_cutouts(kevlar.open(tfile, 'r'))
    target = [record for record in tinstream][0]

    aln = VariantMapping(query, target, 1e6, '25D263M25D')
    variants = aln.call_variants(31)
    assert len(variants) == 2


@pytest.mark.parametrize('targetfile,queryfile,cigar', [
    ('pico-7-refr.fa', 'pico-7-asmbl.fa', '10D83M190D75M20I1M'),
    ('pico-2-refr.fa', 'pico-2-asmbl.fa', '10D89M153I75M20I'),
])
def test_call_formerly_inscrutable(targetfile, queryfile, cigar, capsys):
    target = data_file(targetfile)
    query = data_file(queryfile)
    args = kevlar.cli.parser().parse_args(['call', query, target])
    kevlar.call.main(args)

    out, err = capsys.readouterr()
    print(out)
    assert 'GC=' not in out


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
    tinstream = kevlar.reference.load_refr_cutouts(kevlar.open(tfile, 'r'))
    targetseqs = [record for record in tinstream]

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
    tinstream = kevlar.reference.load_refr_cutouts(kevlar.open(tfile, 'r'))
    target = [record for record in tinstream][0]

    aln = VariantMapping(query, target, 1e6, cigar)
    variants = aln.call_variants(21)
    assert len(variants) == 1
    assert variants[0].window == altwindow
    assert variants[0].refrwindow == refrwindow


def test_nocall():
    # Intentionally mismatched
    qfile = data_file('phony-deletion-01.contig.fa')
    tfile = data_file('phony-insertion-01.gdna.fa')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(qfile, 'r'))
    query = [record for record in qinstream][0]
    tinstream = kevlar.reference.load_refr_cutouts(kevlar.open(tfile, 'r'))
    target = [record for record in tinstream][0]

    aln = VariantMapping(query, target, 1e6, '25D5M22I5M46D8M13D2M35I')
    variants = aln.call_variants(21)
    assert len(variants) == 1
    assert variants[0].vcf == (
        'yourchr\t801\t.\t.\t.\t.\t.\t'
        'CG=25D5M22I5M46D8M13D2M35I;NC=inscrutablecigar;QN=contig4:cc=1;'
        'QS=AACTGGTGGGCTCAAGACTAAAAAGACTTTTTTGGTGACAAGCAGGGCGGCCTGCCCTTCCTGTAG'
        'TGCAAGAAAAT'
    )


@pytest.mark.parametrize('part,coord,window', [
    (12, 7027071, 'CAGGGAGAGGCAGCCTGCCCTCAACCTGGGAGAGCACTGTCTAATCAGCTCCCATCTCA'
                  'GG'),
    (16, 25755121, 'TTTTGGTGTTTAGACATGAAGTCCTTGCCCATCGAGTTATGCCTATGTCCTGAATGCT'
                   'ATTGCCTAGG'),
    (23, 59459928, 'CAGGCGTGAGCCACCGCGCCTGGCCAGGAGCATTGTTTGAACCCAGAAGGCGGAGGTT'
                   'GCA'),
    (192, 28556906, 'AAAATACAAAAATTAGCCAGGCATGGTGGTGCATGCCTGTAATACCAGCCTTTTAGA'
                    'GGC')
])
def test_funky_cigar(part, coord, window):
    contigfile = data_file('funkycigar/part.cc{:d}.contig.fa.gz'.format(part))
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contigfile, 'r'))
    contigs = list(contigstream)

    gdnafile = data_file('funkycigar/part.cc{:d}.gdna.fa.gz'.format(part))
    gdnastream = kevlar.reference.load_refr_cutouts(kevlar.open(gdnafile, 'r'))
    targets = list(gdnastream)

    calls = list(call(targets, contigs))
    assert len(calls) == 1
    assert calls[0].seqid == '17'
    assert calls[0].position == coord - 1
    assert calls[0].info['VW'] == window


def test_funky_cigar_deletion():
    contigfile = data_file('funkycigar/deletion.contig.fa')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contigfile, 'r'))
    contigs = list(contigstream)

    gdnafile = data_file('funkycigar/deletion.gdna.fa')
    gdnastream = kevlar.reference.load_refr_cutouts(kevlar.open(gdnafile, 'r'))
    targets = list(gdnastream)

    calls = list(call(targets, contigs))
    assert len(calls) == 1
    assert calls[0].seqid == 'chr42'
    assert calls[0].position == 53644
    assert calls[0]._refr == 'ATGTCTGTTTTCTTAACCT'
    assert calls[0]._alt == 'A'


def test_perfect_match():
    contigfile = data_file('nodiff.contig.fa')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contigfile, 'r'))
    contigs = list(contigstream)

    gdnafile = data_file('nodiff.gdna.fa')
    gdnastream = kevlar.reference.load_refr_cutouts(kevlar.open(gdnafile, 'r'))
    targets = list(gdnastream)

    calls = list(call(targets, contigs))
    assert len(calls) == 1
    assert calls[0].seqid == 'chr99'
    assert calls[0].position == 2899377
    assert calls[0].info['NC'] == 'perfectmatch'


def test_multibest_revcom():
    contigfile = data_file('multibestrc.contig.fa')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contigfile, 'r'))
    contigs = list(contigstream)

    gdnafile = data_file('multibestrc.gdna.fa')
    gdnastream = kevlar.reference.load_refr_cutouts(kevlar.open(gdnafile, 'r'))
    targets = list(gdnastream)

    calls = list(call(targets, contigs))
    assert len(calls) == 4

    coordinates = [c.position + 1 for c in calls]
    assert coordinates == [34495786, 34583830, 58088279, 60344854]
    for c in calls:
        assert c._refr == 'A'
        assert c._alt == 'G'
        assert c.window == ('CCTGAGCCCTCTCAAGTCGGGTCCTGGCCCGGTCTGCCCATGAGGCTGG'
                            'GCCTGAGCCCCA')


def test_variant_mapping():
    contig = screed.Record(
        name='contig1',
        sequence='CCTGAGCCCTCTCAAGTCGGGTCCTGGCCCGGTCTGCCCATGAGGCTGGGCCTGAGCCCC'
    )
    cutout = kevlar.reference.ReferenceCutout(
        defline='chr1_10000-10060',
        sequence='CCTGAGCCCTCTCAAGTCGGGTCCTGGCCCAGTCTGCCCATGAGGCTGGGCCTGAGCCCC'
    )
    mapping = VariantMapping(contig, cutout, score=1e6, cigar='60M')

    assert mapping.seqid == 'chr1'
    assert mapping.interval == ('chr1', 10000, 10060)


def test_align_mates():
    mate_seqs = kevlar.open(data_file('minitrio/novel-mates.fastq.gz'), 'r')
    record = screed.Record(
        name='bogusread',
        sequence='NNNNN',
        mateseqs=[r.sequence for r in kevlar.parse_augmented_fastx(mate_seqs)]
    )
    refrfile = data_file('minitrio/refr.fa')
    kevlar.reference.autoindex(refrfile)
    positions = list(kevlar.call.align_mates(record, refrfile))
    seqids = set([seqid for seqid, coord in positions])
    coords = sorted([coord for seqid, coord in positions])
    assert seqids == set(['seq1'])
    assert coords == [
        45332, 45377, 45393, 45428, 45440, 45447, 46092, 46093, 46099, 46127,
        46131, 46146, 46148, 48025, 48035
    ]


def test_mate_distance():
    coords = [
        45332, 45377, 45393, 45428, 45440, 45447, 46092, 46093, 46099, 46127,
        46131, 46146, 46148, 48025, 48035
    ]
    positions = [('seq1', c) for c in coords]
    gdna_pos = ('seq1', 45727, 45916)
    assert kevlar.call.mate_distance(positions, gdna_pos) == 506.46666666666664

    positions = [('seq2', 4000), ('seq2', 3000), ('seq2', 5100), ('seq3', 1)]
    gdna_pos = ('seq2', 5000, 5500)
    assert kevlar.call.mate_distance(positions, gdna_pos) == 1000.0


@pytest.mark.parametrize('query,target,dist,n,msgcount', [
    ('phony-snv-01b.contig.fa', 'phony-snv-01.gdna.fa', 5, 1, 1),
    ('phony-snv-02b.contig.fa', 'phony-snv-02.gdna.fa', 5, 1, 1),
    ('phony-snv-01b.contig.fa', 'phony-snv-01.gdna.fa', 2, 2, 0),
    ('phony-snv-02b.contig.fa', 'phony-snv-02.gdna.fa', None, 2, 0),
])
def test_call_near_end(query, target, dist, n, msgcount):
    log = StringIO()
    contig = next(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file(query), 'r')
        )
    )
    cutout = next(
        kevlar.reference.load_refr_cutouts(
            kevlar.open(data_file(target), 'r')
        )
    )
    aln = kevlar.call.align_both_strands(cutout, contig)
    calls = list(aln.call_variants(31, mindist=dist, logstream=log))
    assert len(calls) == n

    err = log.getvalue()
    count = err.count('discarding SNV due to proximity to end of the contig')
    assert count == msgcount


def test_call_truncated_windows_snv():
    contig = next(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file('trunc-snv.contig.fa'), 'r')
        )
    )
    cutout = next(
        kevlar.reference.load_refr_cutouts(
            kevlar.open(data_file('trunc-snv.gdna.fa'), 'r')
        )
    )
    aln = kevlar.call.align_both_strands(cutout, contig)
    calls = list(aln.call_variants(31))
    assert len(calls) == 1
    assert calls[0].window == 'TAGCATACAGGTAGTCAGGGGGTGTCTGCGACCACAGCTGAA'
    assert calls[0].refrwindow == 'TAGCATACAGGAAGTCAGGGGGTGTCTGCGACCACAGCTGAA'


def test_call_truncated_windows_indel():
    contig = next(
        kevlar.parse_augmented_fastx(
            kevlar.open(data_file('trunc-indel.contig.augfasta'), 'r')
        )
    )
    cutout = next(
        kevlar.reference.load_refr_cutouts(
            kevlar.open(data_file('trunc-indel.gdna.augfasta'), 'r')
        )
    )
    aln = kevlar.call.align_both_strands(cutout, contig)
    calls = list(aln.call_variants(31))
    assert len(calls) == 1
    assert calls[0].window == 'GTATATACACATATATGTGTATATATGTGTATATATGT'
    assert calls[0].refrwindow == ('GTATATACATATATGTGTATATATGTGTATATACATATATGT'
                                   'GTATATATGTGTATATATGT')
