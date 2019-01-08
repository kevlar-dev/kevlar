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
from . import count
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
from . import varfilter
from . import simlike
from . import split
from . import dist
from . import unband

mains = {
    'count': kevlar.count.main,
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
    'varfilter': kevlar.varfilter.main,
    'simlike': kevlar.simlike.main,
    'split': kevlar.split.main,
    'dist': kevlar.dist.main,
    'unband': kevlar.unband.main,
}

subparser_funcs = {
    'count': count.subparser,
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
    'varfilter': varfilter.subparser,
    'simlike': simlike.subparser,
    'split': split.subparser,
    'dist': dist.subparser,
    'unband': unband.subparser,
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
    parser.add_argument('-l', '--logfile', metavar='F', help='log file for '
                        'diagnostic messages, warnings, and errors')
    parser.add_argument('--tee', action='store_true', help='write diagnostic '
                        'output to logfile AND terminal (stderr)')
    subparsers = parser.add_subparsers(dest='cmd', metavar='cmd',
                                       help='"' + subcommandstr + '"')
    for func in subparser_funcs.values():
        func(subparsers)

    return parser


def parse_args(arglist=None):
    args = kevlar.cli.parser().parse_args(arglist)
    kevlar.logstream = sys.stderr
    if args.logfile and args.logfile != '-':
        kevlar.logstream = kevlar.open(args.logfile, 'w')
    kevlar.teelog = args.tee
    return args
