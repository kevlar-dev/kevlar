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
import shutil
import tempfile
import kevlar
import kevlar.__main__


@pytest.mark.long
def test_trio2(capsys):
    from sys import stdout, stderr
    tempdir = tempfile.mkdtemp()
    novelouts = ['{}/out{}'.format(tempdir, i) for i in range(4)]
    case = kevlar.tests.data_file('trio2/case1.fq.gz')
    controls = kevlar.tests.data_glob('trio2/ctrl[1,2].fq.gz')
    for i in range(4):
        arglist = [
            'novel', '--case', case,
            '--control', controls[0], '--control', controls[1],
            '--band', str(i+1), '--num-bands', '4', '--out', novelouts[i],
            '--memory', '200K', '--ksize', '31',
            '--case-min', '8', '--ctrl-max', '1'
        ]
        args = kevlar.cli.parser().parse_args(arglist)
        kevlar.novel.main(args)

    arglist = [
        'collect', '--memory', '5K', '--ksize', '31', '--minabund', '8',
        '--collapse'
    ] + novelouts
    args = kevlar.cli.parser().parse_args(arglist)
    args.out = stdout
    args.logfile = stderr
    kevlar.__main__.main(args)
    out, err = capsys.readouterr()

    assert '1 collapsed linear paths' in err
    lastline = out.strip().split('\n')[-1]
    contig, contigrc = lastline.split(',')[0:2]
    assert 'AGCCTCTG' in contig or 'AGCCTCTG' in contigrc

    shutil.rmtree(tempdir)
