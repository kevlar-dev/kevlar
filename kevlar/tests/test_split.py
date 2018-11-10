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
from kevlar.tests import data_file
from tempfile import mkdtemp
from shutil import rmtree


def test_split():
    infile = data_file('fiveparts.augfastq.gz')
    partfile = kevlar.open(infile, 'r')
    readstream = kevlar.parse_augmented_fastx(partfile)
    partstream = kevlar.parse_partitioned_reads(readstream)
    outstreams = [StringIO(), StringIO(), StringIO()]
    kevlar.split.split(partstream, outstreams)

    assert 'kvcc=1' in outstreams[0].getvalue()
    assert 'kvcc=2' in outstreams[1].getvalue()
    assert 'kvcc=3' in outstreams[2].getvalue()
    assert 'kvcc=4' in outstreams[0].getvalue()
    assert 'kvcc=5' in outstreams[1].getvalue()


def test_split_cli():
    infile = data_file('fiveparts.augfastq.gz')
    tempdir = mkdtemp()
    print(tempdir)
    arglist = ['split', infile, '3', tempdir + '/out']
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.split.main(args)

    outfile = tempdir + '/out.1'
    readstream = kevlar.parse_augmented_fastx(kevlar.open(outfile, 'r'))
    partstream = kevlar.parse_partitioned_reads(readstream)
    partitions = list(partstream)
    partitions = [part for partid, part in partitions]
    assert len(partitions) == 2
    assert len(partitions[0]) == 67
    assert len(partitions[1]) == 12

    outfile = tempdir + '/out.2'
    readstream = kevlar.parse_augmented_fastx(kevlar.open(outfile, 'r'))
    partstream = kevlar.parse_partitioned_reads(readstream)
    partitions = list(partstream)
    partitions = [part for partid, part in partitions]
    assert len(partitions) == 2
    assert len(partitions[0]) == 23
    assert len(partitions[1]) == 11

    outfile = tempdir + '/out.3'
    readstream = kevlar.parse_augmented_fastx(kevlar.open(outfile, 'r'))
    partstream = kevlar.parse_partitioned_reads(readstream)
    partitions = list(partstream)
    partitions = [part for partid, part in partitions]
    assert len(partitions) == 1
    assert len(partitions[0]) == 15

    rmtree(tempdir)
