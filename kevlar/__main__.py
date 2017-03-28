#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import argparse
import kevlar


def main(args=None):
    """
    Entry point for the tag CLI.

    Isolated as a method so that the CLI can be called by other Python code
    (e.g. for testing), in which case the arguments are passed to the function.
    If no arguments are passed to the function, parse them from the command
    line.
    """
    if args is None:
        args = kevlar.cli.parser().parse_args()

    if args.cmd is None:
        kevlar.cli.parser().parse_args(['-h'])

    assert args.cmd in kevlar.cli.mains
    mainmethod = mains[args.cmd]
    mainmethod(args)
