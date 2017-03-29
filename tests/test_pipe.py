#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import re
import shutil
import tempfile
import kevlar


@pytest.mark.long
def test_trio2(capsys):
    from sys import stdout, stderr
    tempdir = tempfile.mkdtemp()
    findouts = ['{}/out{}'.format(tempdir, i) for i in range(4)]
    for i in range(4):
        arglist = [
            'find', '--cases', 'tests/data/trio2/case1.fq.gz',
            '--band', str(i+1), '--num-bands', '4', '--out', findouts[i],
            '--memory', '200K', '--ksize', '31', '--case_min', '8',
            '--ctrl_max', '1', '--controls', 'tests/data/trio2/ctrl1.fq.gz',
            'tests/data/trio2/ctrl2.fq.gz'
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.find.main(args)

    arglist = [
        'collect', '--memory', '5K', '--ksize', '31', '--minabund', '8',
        '--collapse'
    ] + findouts
    args = kevlar.cli.parser().parse_args(arglist)
    args.out = stdout
    args.logfile = stderr
    kevlar.collect.main(args)
    out, err = capsys.readouterr()

    assert '1 collapsed linear paths' in err
    lastline = out.strip().split('\n')[-1]
    contig, contigrc = lastline.split(',')[0:2]
    assert 'AGCCTCTG' in contig or 'AGCCTCTG' in contigrc

    shutil.rmtree(tempdir)
