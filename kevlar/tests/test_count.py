#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import glob
import os
import pytest
import re
from tempfile import NamedTemporaryFile
import screed
import time
import kevlar
from khmer import Nodetable, _buckets_per_byte
from kevlar.tests import data_file, data_glob


@pytest.fixture
def triomask():
    mask = khmer.Counttable(19, 1e4, 4)
    mask.consume('TGAGGGGACTAGGTGATCAGGTGAGGGTTTCCCAGTTCCCGAAGATGACT')
    mask.consume('GATCTTTCGCTCCCTGTCATCAAGGAGTGATACGCGAAGTGCGTCCCCTT')
    mask.consume('GAAGTTTTGACAATTTACGTGAGCCCTACCTAGCGAAACAACAGAGAACC')
    return mask


@pytest.mark.parametrize('mask,numbands,band', [
    (None, None, None),
    (None, 9, 2),
    (triomask, None, None),
    (triomask, 23, 19),
])
def test_load_threading(mask, numbands, band):
    # Smoke test: make sure things don't explode when run in "threaded" mode.
    infiles = data_glob('trio1/case1.fq')
    sketch = kevlar.count.load_sample_seqfile(
        infiles, 19, 1e7, mask=mask, numbands=numbands, band=band, numthreads=2
    )


@pytest.mark.parametrize('infile,testout,numbands,band,kmers_stored', [
    ('case', 'case', 0, 0, 973),
    ('ctrl1', 'ctrl1', 0, 0, 973),
    ('ctrl2', 'ctrl2', 0, 0, 966),
    ('case', 'case-band-2-1', 2, 1, 501),
    ('case', 'case-band-16-7', 16, 7, 68),
])
def test_count_simple(infile, testout, numbands, band, kmers_stored, capsys):
    infile = data_file('simple-genome-{}-reads.fa.gz'.format(infile))
    testout = data_file('simple-genome-{}.ct'.format(testout))
    with NamedTemporaryFile() as outfile:
        arglist = ['count', '--ksize', '25', '--memory', '10K',
                   '--num-bands', str(numbands), '--band', str(band),
                   outfile.name, infile]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.count.main(args)
        out, err = capsys.readouterr()

        assert '600 reads processed' in str(err)
        assert '{:d} distinct k-mers stored'.format(kmers_stored) in str(err)

        outputfilename = outfile.name + '.counttable'
        with open(outputfilename, 'rb') as f1, open(testout, 'rb') as f2:
            assert f1.read() == f2.read()


def test_count_threading():
    with NamedTemporaryFile(suffix='.counttable') as outfile:
        infile = data_file('trio1/case1.fq')
        arglist = ['count', '--ksize', '19', '--memory', '500K',
                   '--threads', '2', outfile.name, infile]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.count.main(args)

    # No checks, just doing a "smoke test" to make sure things don't explode
    # when counting is done in "threaded" mode.


def test_count_problematic():
    arglist = [
        'count', '--ksize', '21', '--memory', '200K', '--band', '2',
        'bogusoutput', data_file('trio1/ctrl1.fq')
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    with pytest.raises(ValueError) as ve:
        kevlar.count.main(args)
    assert 'Must specify --num-bands and --band together' in str(ve)

    arglist = [
        'count', '--ksize', '21', '--memory', '97',
        'bogusoutput', data_file('trio1/ctrl1.fq')
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    with pytest.raises(kevlar.sketch.KevlarUnsuitableFPRError):
        kevlar.count.main(args)


def test_effcount_smoketest():
    with NamedTemporaryFile(suffix='.ct') as o1, \
            NamedTemporaryFile(suffix='.ct') as o2, \
            NamedTemporaryFile(suffix='.ct') as o3:

        arglist = [
            'effcount', '--sample', data_file('trio1/ctrl1.fq'),
            '--sample', data_file('trio1/ctrl2.fq'),
            '--sample', data_file('trio1/case2.fq'),
            '--ksize', '21', '--memory', '200K', '--memfrac', '0.005',
            '--max-fpr', '0.1', '--max-abund', '0', '--threads', '2',
            o1.name, o2.name, o3.name
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.effcount.main(args)


def test_effcount_problematic():
    arglist = [
        'effcount', '--sample', data_file('trio1/ctrl1.fq'), '--ksize', '21',
        '--memory', '200K', 'bogusoutput'
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    with pytest.raises(AssertionError):
        kevlar.effcount.main(args)

    arglist = [
        'effcount', '--sample', data_file('trio1/ctrl1.fq'),
        '--sample', data_file('trio1/ctrl2.fq'), '--ksize', '21',
        '--memory', '200K', 'bogusoutput'
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    with pytest.raises(ValueError) as ve:
        kevlar.effcount.main(args)
    assert 'number of outfiles must match number of declared' in str(ve)

    arglist = [
        'effcount', '--sample', data_file('trio1/ctrl1.fq'),
        '--sample', data_file('trio1/ctrl2.fq'), '--ksize', '21',
        '--memory', '200K', '--band', '2', 'bogusoutput1', 'bogusoutput2'
    ]
    args = kevlar.cli.parser().parse_args(arglist)
    with pytest.raises(ValueError) as ve:
        kevlar.effcount.main(args)
    assert 'Must specify --num-bands and --band together' in str(ve)


@pytest.mark.parametrize('count,smallcount,extension,shortext', [
    (True, True, '.smallcounttable', '.sct'),
    (True, False, '.counttable', '.ct'),
    (False, True, '.nodetable', '.nt'),
    (False, False, '.nodetable', '.nt'),
])
def test_load_sample_seqfile(count, smallcount, extension, shortext):
    infile = data_file('bogus-genome/refr.fa')
    with NamedTemporaryFile() as outfile:
        sketch = kevlar.count.load_sample_seqfile(
            [infile], 21, 1e6, count=count, smallcount=smallcount,
            outfile=outfile.name
        )
        assert sketch.get('GAATCGGTGGCTGGTTGCCGT') > 0
        assert sketch.get('GATTACAGATTACAGATTACA') == 0
        assert os.path.exists(outfile.name + extension)

    with NamedTemporaryFile(suffix=shortext) as outfile:
        sketch = kevlar.count.load_sample_seqfile(
            [infile], 21, 1e6, count=count, smallcount=smallcount,
            outfile=outfile.name
        )
        assert sketch.get('GAATCGGTGGCTGGTTGCCGT') > 0
        assert sketch.get('GATTACAGATTACAGATTACA') == 0
        assert not os.path.exists(outfile.name + extension)
        assert os.path.exists(outfile.name)


@pytest.mark.parametrize('count,smallcount,count_masked,kpresent,kabsent', [
    (True, True, True, 'CACCAATCCGTACGGAGAGCC', 'GAATCGGTGGCTGGTTGCCGT'),
    (True, False, True, 'CACCAATCCGTACGGAGAGCC', 'GAATCGGTGGCTGGTTGCCGT'),
    (False, True, True, 'CACCAATCCGTACGGAGAGCC', 'GAATCGGTGGCTGGTTGCCGT'),
    (False, False, True, 'CACCAATCCGTACGGAGAGCC', 'GAATCGGTGGCTGGTTGCCGT'),
    (True, True, False, 'GAATCGGTGGCTGGTTGCCGT', 'CACCAATCCGTACGGAGAGCC'),
    (True, False, False, 'GAATCGGTGGCTGGTTGCCGT', 'CACCAATCCGTACGGAGAGCC'),
    (False, True, False, 'GAATCGGTGGCTGGTTGCCGT', 'CACCAATCCGTACGGAGAGCC'),
    (False, False, False, 'GAATCGGTGGCTGGTTGCCGT', 'CACCAATCCGTACGGAGAGCC'),
])
def test_load_sample_seqfile_withmask(count, smallcount, count_masked,
                                      kpresent, kabsent):
    mask = Nodetable(21, 1e4, 4)
    mask.consume('CACCAATCCGTACGGAGAGCCGTATATATAGACTGCTATACTATTGGATCGTACGGGGC')
    sketch = kevlar.count.load_sample_seqfile(
        [data_file('bogus-genome/refr.fa')], 21, 1e6, mask=mask,
        consume_masked=count_masked, count=count, smallcount=smallcount,
    )
    assert sketch.get(kpresent) > 0
    assert sketch.get(kabsent) == 0
    assert sketch.get('GATTACAGATTACAGATTACA') == 0


def test_count_cli_with_mask(capsys):
    mask = Nodetable(21, 1e4, 4)
    mask.consume('CACCAATCCGTACGGAGAGCCGTATATATAGACTGCTATACTATTGGATCGTACGGGGC')
    with NamedTemporaryFile(suffix='.nt') as maskfile, \
            NamedTemporaryFile(suffix='.sct') as countfile:
        mask.save(maskfile.name)
        arglist = [
            'count', '--ksize', '21', '--mask', maskfile.name,
            '--memory', '1M', countfile.name, data_file('bogus-genome/refr.fa')
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.count.main(args)
    out, err = capsys.readouterr()
    assert '36898 distinct k-mers stored' in err


@pytest.mark.parametrize('count,smallcount,sketchtype', [
    (False, False, 'nodegraph'),
    (True, False, 'countgraph'),
    (True, True, 'smallcountgraph'),
])
def test_load_sample_seqfile_memory_test(count, smallcount, sketchtype):
    requested_memory = 2e6
    sketch = kevlar.count.load_sample_seqfile(
        [data_file('bogus-genome/refr.fa')], 21, requested_memory, count=count,
        smallcount=smallcount,
    )
    buckets = sum(sketch.hashsizes())
    actual_memory = buckets / _buckets_per_byte[sketchtype]
    assert actual_memory / requested_memory == pytest.approx(1.0, rel=1e-4)
