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
