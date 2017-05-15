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
import kevlar
import sys


mains = {
    'dump': kevlar.dump.main,
    'count': kevlar.count.main,
    'novel': kevlar.novel.main,
    'collect': kevlar.collect.main,
    'filter': kevlar.filter.main,
    'reaugment': kevlar.reaugment.main,
    'assemble': kevlar.assemble.main,
    'mutate': kevlar.mutate.main,
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
    kevlar.dump.subparser(subparsers)
    kevlar.count.subparser(subparsers)
    kevlar.novel.subparser(subparsers)
    kevlar.collect.subparser(subparsers)
    kevlar.filter.subparser(subparsers)
    kevlar.reaugment.subparser(subparsers)
    kevlar.assemble.subparser(subparsers)
    kevlar.mutate.subparser(subparsers)

    return parser
