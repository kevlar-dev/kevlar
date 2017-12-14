#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------


class MutableString(object):
    """
    Mutable string class

    Borrows heavily from https://stackoverflow.com/a/10572792/459780.
    """

    def __init__(self, data):
        self.data = list(data)

    def __str__(self):
        return ''.join(self.data)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return str(self) == str(other)

    def __add__(self, chars):
        """
        Addition operator

        For consistency with other Python objects, the addition operator
        creates a new object rather than appending in place. However, this
        makes a copy of the data which goes against the performance this data
        structure optimizes for.
        """
        newdata = list(self.data) + list(str(chars))
        return MutableString(''.join(newdata))

    def __iadd__(self, chars):
        self.data.extend(list(str(chars)))
        return self

    def __contains__(self, teststr):
        return teststr in str(self)

    def __setitem__(self, index, value):
        self.data[index] = value

    def __getitem__(self, index):
        if type(index) == slice:
            return ''.join(self.data[index])
        return self.data[index]

    def __delitem__(self, index):
        del self.data[index]

    def __len__(self):
        return len(self.data)
