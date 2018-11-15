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
from . import effcount
from . import novel
from . import filter
from . import augment
from . import assemble
from . import mutate
from . import gentrio
from . import partition
from . import localize
from . import call
from . import alac
from . import simplex
from . import simlike
from . import split
from . import dist

mains = {
    'dump': kevlar.dump.main,
    'count': kevlar.count.main,
    'effcount': kevlar.effcount.main,
    'novel': kevlar.novel.main,
    'filter': kevlar.filter.main,
    'augment': kevlar.augment.main,
    'assemble': kevlar.assemble.main,
    'mutate': kevlar.mutate.main,
    'gentrio': kevlar.gentrio.main,
    'partition': kevlar.partition.main,
    'localize': kevlar.localize.main,
    'call': kevlar.call.main,
    'alac': kevlar.alac.main,
    'simplex': kevlar.simplex.main,
    'simlike': kevlar.simlike.main,
    'split': kevlar.split.main,
    'dist': kevlar.dist.main,
}

subparser_funcs = {
    'dump': dump.subparser,
    'count': count.subparser,
    'effcount': effcount.subparser,
    'novel': novel.subparser,
    'filter': filter.subparser,
    'augment': augment.subparser,
    'assemble': assemble.subparser,
    'mutate': mutate.subparser,
    'gentrio': gentrio.subparser,
    'partition': partition.subparser,
    'localize': localize.subparser,
    'call': call.subparser,
    'alac': alac.subparser,
    'simplex': simplex.subparser,
    'simlike': simlike.subparser,
    'split': split.subparser,
    'dist': dist.subparser,
}


def parser():
    bubbletext = r'''
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
┌ ┐            ┌ ┐
| |            | |
| | _______   _| | __ _ _ __
| |/ / _ \ \ / / |/ _` | '__|
|   <  __/\ V /| | (_| | |        reference-free variant discovery in
|_|\_\___| \_/ |_|\__,_|_|                   large eukaryotic genomes
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
'''
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
