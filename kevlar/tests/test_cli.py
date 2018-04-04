#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import kevlar


def test_kevlar_open():
    thefile = kevlar.tests.data_file('wasp-pass.contig.augfasta')
    filehandle = kevlar.open(thefile, 'r')
    filecontents = filehandle.read()
    assert len(filecontents.strip().split('\n')) == 9

    with pytest.raises(ValueError) as ve:
        filehandle = kevlar.open(thefile, 'p')
    assert 'invalid mode "p"' in str(ve)


def test_main(capsys):
    import kevlar.__main__
    contig = kevlar.tests.data_file('wasp-pass.contig.augfasta')
    cutout = kevlar.tests.data_file('wasp.gdna.fa')
    arglist = ['call', contig, cutout]
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.__main__.main(args)
    out, err = capsys.readouterr()
    assert '\tPASS\t' in out
    assert '\tPassengerVariant\t' in out


def test_help(capsys):
    with pytest.raises(SystemExit):
        kevlar.cli.parser().parse_args(['-h'])
    out, err = capsys.readouterr()
    assert 'show this help message and exit' in out


def test_version(capsys):
    with pytest.raises(SystemExit):
        kevlar.cli.parser().parse_args(['-v'])
    out, err = capsys.readouterr()
    assert kevlar.__version__ in out or kevlar.__version__ in err


@pytest.mark.parametrize('subcommand', [s for s in kevlar.cli.mains])
def test_help_sub(subcommand, capsys):
    with pytest.raises(SystemExit):
        kevlar.cli.parser().parse_args([subcommand, '-h'])
    out, err = capsys.readouterr()
    assert subcommand in out
    assert 'show this help message and exit' in out
