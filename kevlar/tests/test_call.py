#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


import filecmp
import kevlar
from kevlar.call import call
from kevlar.sequence import Record
from kevlar.tests import data_file
import pytest
from tempfile import NamedTemporaryFile


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
    """Ensure isolated SNVs are called correctly.

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
    print('DEBUG', calls[0].vcf)
    assert calls[0].seqid == '17'
    assert calls[0].position == coord - 1
    assert calls[0].attribute('ALTWINDOW') == window


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
    assert calls[0].attribute('CONTIG') == contigs[0].sequence


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
    assert calls[0].filterstr == 'PerfectMatch'


def test_cigar_filter_regression():
    contigfile = data_file('14153.cc5463.contig.augfasta.gz')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contigfile, 'r'))
    contigs = list(contigstream)

    gdnafile = data_file('14153.cc5463.gdna.augfasta.gz')
    gdnastream = kevlar.reference.load_refr_cutouts(kevlar.open(gdnafile, 'r'))
    targets = list(gdnastream)

    calls = sorted(call(targets, contigs), key=lambda c: c.position)
    assert len(calls) == 2
    assert calls[1].seqid == '6'

    # Equally valid calls from equally optimal alignments
    c1 = ('AGAAA', 'A', 154734241)
    c2 = ('GAAGA', 'G', 154734239)

    varcall = (calls[1]._refr, calls[1]._alt, calls[1].position)
    assert varcall in (c1, c2)


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


def test_snv_dedup():
    contigfile = data_file('bee-dupl.contigs.augfasta')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contigfile, 'r'))
    contigs = list(contigstream)

    gdnafile = data_file('bee-dupl.gdna.fa')
    gdnastream = kevlar.reference.load_refr_cutouts(kevlar.open(gdnafile, 'r'))
    targets = list(gdnastream)

    calls = list(call(targets, contigs, ksize=27))
    assert len(calls) == 1
    assert calls[0].seqid == 'linkagegroup5'
    assert calls[0].position == 8174 - 1


def test_debug_mode(capsys):
    contig = data_file('wasp-pass.contig.augfasta')
    gdna = data_file('wasp.gdna.fa')
    arglist = ['call', '--debug', contig, gdna]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.call.main(args)

    out, err = capsys.readouterr()
    alignstr = kevlar.open(data_file('wasp-align.txt'), 'r').read().strip()
    assert alignstr in err


def test_call_generate_mask():
    contigfile = data_file('fiveparts.contigs.augfasta.gz')
    gdnafile = data_file('fiveparts.gdnas.fa.gz')
    refrfile = data_file('fiveparts-refr.fa.gz')
    with NamedTemporaryFile(suffix='.nt') as maskfile:
        arglist = [
            'call', '--gen-mask', maskfile.name, '--mask-mem', '1M',
            '--refr', refrfile, contigfile, gdnafile
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.call.main(args)
        testfilename = data_file('fiveparts-genmask.nodetable')
        assert filecmp.cmp(testfilename, maskfile.name) is True


def test_call_generate_mask_lowmem(capsys):
    contigfile = data_file('fiveparts.contigs.augfasta.gz')
    gdnafile = data_file('fiveparts.gdnas.fa.gz')
    refrfile = data_file('fiveparts-refr.fa.gz')
    with NamedTemporaryFile(suffix='.nt') as maskfile:
        arglist = [
            'call', '--gen-mask', maskfile.name, '--mask-mem', '100',
            '--refr', refrfile, contigfile, gdnafile
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.call.main(args)
    out, err = capsys.readouterr()
    message = 'WARNING: mask FPR is 0.8065; exceeds user-specified limit'
    assert message in out or message in err


def test_call_mnv():
    contigfile = data_file('mnv-contig.augfasta')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contigfile, 'r'))
    contigs = list(contigstream)

    gdnafile = data_file('mnv-gdna.fa')
    gdnastream = kevlar.reference.load_refr_cutouts(kevlar.open(gdnafile, 'r'))
    targets = list(gdnastream)

    caller = kevlar.call.call(targets, contigs, ksize=49)
    calls = list(caller)
    calls.sort(key=lambda v: v.position)

    assert len(calls) == 3
    assert [v.position for v in calls] == [98153308, 98153312, 98153407]
    assert calls[1]._refr == 'GA'
    assert calls[1]._alt == 'TT'
    assert calls[2].filterstr == 'PassengerVariant'


def test_call_mnv_3bp():
    contigfile = data_file('ant.contig.augfasta')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contigfile, 'r'))
    contigs = list(contigstream)

    gdnafile = data_file('ant.gdna.fa')
    gdnastream = kevlar.reference.load_refr_cutouts(kevlar.open(gdnafile, 'r'))
    targets = list(gdnastream)

    caller = kevlar.call.call(targets, contigs, ksize=29)
    calls = list(caller)

    assert len(calls) == 1
    assert calls[0]._refr == 'ACG'
    assert calls[0]._alt == 'GTT'
    assert calls[0].filterstr == 'PASS'


def test_call_homopolymers():
    contigfile = data_file('homopolymer/14153-6parts.contigs.augfasta')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contigfile, 'r'))
    contigs = list(contigstream)

    gdnafile = data_file('homopolymer/14153-6parts.targets.fasta')
    gdnastream = kevlar.reference.load_refr_cutouts(kevlar.open(gdnafile, 'r'))
    targets = list(gdnastream)

    caller = kevlar.call.call(targets, contigs, ksize=49)
    calls = list(caller)

    assert len(calls) == 6
    filters = [c.filterstr for c in calls]
    assert 'PASS' not in filters
    for f in filters:
        assert 'Homopolymer' in f


def test_call_homopolymers_mixed_results():
    contigfile = data_file('homopolymer/12175-3parts.contigs.augfasta')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contigfile, 'r'))
    partstream = kevlar.parse_partitioned_reads(contigstream)
    contigs = kevlar.call.load_contigs(partstream)

    gdnafile = data_file('homopolymer/12175-3parts.targets.fasta')
    gdnastream = kevlar.reference.load_refr_cutouts(kevlar.open(gdnafile, 'r'))
    partstream = kevlar.parse_partitioned_reads(gdnastream)
    targets = kevlar.call.load_contigs(partstream)

    prelimcalls = list()
    for partid in contigs:
        contiglist = contigs[partid]
        gdnalist = targets[partid]
        caller = kevlar.call.call(gdnalist, contiglist, partid=partid)
        prelimcalls.extend(list(caller))

    kid = kevlar.sketch.load(data_file('homopolymer/12175-kid.sct'))
    mom = kevlar.sketch.load(data_file('homopolymer/12175-mom.sct'))
    dad = kevlar.sketch.load(data_file('homopolymer/12175-dad.sct'))
    refr = kevlar.sketch.load(data_file('homopolymer/12175-refr.sct'))
    scorer = kevlar.simlike.simlike(
        prelimcalls, kid, [mom, dad], refr,
        samplelabels=['Proband', 'Mother', 'Father'],
    )
    calls = list(scorer)

    assert len(calls) == 6
    for c in calls:
        print(c.vcf)
    unintrstng = [c for c in calls if c.filterstr in ('PASS', 'Homopolymer')]
    assert len(unintrstng) == 3

    call1, call2, call3 = unintrstng
    assert call1.position == 123651924
    assert call1.filterstr == 'PASS'  # negative control
    assert call1._refr == 'TAA'
    assert call1._alt == 'T'
    assert call2.position == 124641259
    assert call2.filterstr == 'PASS'  # borderline
    assert call2._refr == 'TAAA'
    assert call2._alt == 'T'
    assert call3.position == 128660727
    assert call3.filterstr == 'Homopolymer'  # positive control


def test_call_homopolymer_filter_disabled():
    contigfile = data_file('homopolymer/12175-3parts.contigs.augfasta')
    contigstream = kevlar.parse_augmented_fastx(kevlar.open(contigfile, 'r'))
    partstream = kevlar.parse_partitioned_reads(contigstream)
    contigs = kevlar.call.load_contigs(partstream)

    gdnafile = data_file('homopolymer/12175-3parts.targets.fasta')
    gdnastream = kevlar.reference.load_refr_cutouts(kevlar.open(gdnafile, 'r'))
    partstream = kevlar.parse_partitioned_reads(gdnastream)
    targets = kevlar.call.load_contigs(partstream)

    prelimcalls = list()
    for partid in contigs:
        contiglist = contigs[partid]
        gdnalist = targets[partid]
        caller = kevlar.call.call(
            gdnalist, contiglist, partid=partid, homopolyfilt=False
        )
        prelimcalls.extend(list(caller))

    kid = kevlar.sketch.load(data_file('homopolymer/12175-kid.sct'))
    mom = kevlar.sketch.load(data_file('homopolymer/12175-mom.sct'))
    dad = kevlar.sketch.load(data_file('homopolymer/12175-dad.sct'))
    refr = kevlar.sketch.load(data_file('homopolymer/12175-refr.sct'))
    scorer = kevlar.simlike.simlike(
        prelimcalls, kid, [mom, dad], refr,
        samplelabels=['Proband', 'Mother', 'Father'],
    )
    calls = list(scorer)

    assert len(calls) == 6
    for c in calls:
        assert 'Homopolymer' not in c.filterstr
