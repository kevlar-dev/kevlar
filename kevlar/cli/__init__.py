#!/usr/bin/env python
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
    'find': kevlar.find.main,
    'collect': kevlar.collect.main,
    'filter': kevlar.filter.main,
}


def parser():
    subcommandstr = '", "'.join(sorted(list(mains.keys())))
    parser = argparse.ArgumentParser()
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
    kevlar.find.subparser(subparsers)
    kevlar.collect.subparser(subparsers)
    kevlar.filter.subparser(subparsers)

    return parser
