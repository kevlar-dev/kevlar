#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import glob
import pkg_resources


def data_file(basename):
    datadir = pkg_resources.resource_filename('kevlar', 'tests/data')
    return datadir + '/' + basename


def data_glob(globstr):
    datadir = pkg_resources.resource_filename('kevlar', 'tests/data')
    return sorted(glob.glob(datadir + '/' + globstr))
