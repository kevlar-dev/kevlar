#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
import sys
import kevlar
from . import dump
from . import count
from . import novel
from . import collect
from . import filter
from . import reaugment
from . import assemble
from . import mutate
from . import partition
from . import localize
from . import call

mains = {
    'dump': kevlar.dump.main,
    'count': kevlar.count.main,
    'novel': kevlar.novel.main,
    'collect': kevlar.collect.main,
    'filter': kevlar.filter.main,
    'reaugment': kevlar.reaugment.main,
    'assemble': kevlar.assemble.main,
    'mutate': kevlar.mutate.main,
    'partition': kevlar.partition.main,
    'localize': kevlar.localize.main,
    'call': kevlar.call.main,
}

subparser_funcs = {
    'dump': dump.subparser,
    'count': count.subparser,
    'novel': novel.subparser,
    'collect': collect.subparser,
    'filter': filter.subparser,
    'reaugment': reaugment.subparser,
    'assemble': assemble.subparser,
    'mutate': mutate.subparser,
    'partition': partition.subparser,
    'localize': localize.subparser,
    'call': call.subparser,
}


def parser():
    bubbletext = """
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
┌ ┐            ┌ ┐
| |            | |
| | _______   _| | __ _ _ __
| |/ / _ \ \ / / |/ _` | '__|
|   <  __/\ V /| | (_| | |        reference-free variant discovery in
|_|\_\___| \_/ |_|\__,_|_|                   large eukaryotic genomes
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
"""
    subcommandstr = '", "'.join(sorted(list(mains.keys())))
    parser = argparse.ArgumentParser(
        description=bubbletext,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser._positionals.title = 'Subcommands'
    parser._optionals.title = 'Global arguments'
    parser.add_argument('-v', '--version', action='version',
                        version='kevlar v{}'.format(kevlar.__version__))
    parser.add_argument('-l', '--logfile', metavar='F', default=sys.stderr,
                        type=argparse.FileType('w'), help='log file for '
                        'diagnostic messages, warnings, and errors')
    subparsers = parser.add_subparsers(dest='cmd', metavar='cmd',
                                       help='"' + subcommandstr + '"')
    for func in subparser_funcs.values():
        func(subparsers)

    return parser
