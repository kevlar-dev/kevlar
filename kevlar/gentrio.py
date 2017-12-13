#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict, namedtuple
import random
import sys
import kevlar


nucl_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
index_to_nucl = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}


def mutate_snv(sequence, position, offset, ksize=31):
    orignucl = sequence[position]
    nuclindex = nucl_to_index[orignucl]
    newindex = (nuclindex + offset) % 4
    newnucl = index_to_nucl[newindex]

    windowstart = max(position - ksize + 1, 0)
    windowend = min(position + ksize, len(sequence))
    refrwindow = sequence[windowstart:windowend]
    altwindow = '{:s}{:s}{:s}'.format(
        sequence[windowstart:position], newnucl,
        sequence[position + 1:windowend]
    )

    return orignucl, newnucl, refrwindow, altwindow


def mutate_insertion(sequence, position, length, duplpos, ksize=31):
    duplseq = sequence[duplpos:duplpos + length]
    refrseq = sequence[position - 1]
    altseq = refrseq + duplseq

    windowstart = max(position - ksize + 1, 0)
    windowend = min(position + ksize - 1, len(sequence))
    refrwindow = sequence[windowstart:windowend]
    altwindow = '{:s}{:s}{:s}'.format(
        sequence[windowstart:position], duplseq, sequence[position:windowend]
    )

    return refrseq, altseq, refrwindow, altwindow


def mutate_deletion(sequence, position, length, ksize=31):
    delseq = sequence[position:position + length]
    altseq = sequence[position - 1]
    refrseq = altseq + delseq

    windowstart = max(position - ksize + 1, 0)
    windowend = min(position + length + ksize - 1, len(sequence))
    refrwindow = sequence[windowstart:windowend]
    altwindow = '{:s}{:s}'.format(
        sequence[windowstart:position], sequence[position + length:windowend]
    )

    return refrseq, altseq, refrwindow, altwindow


def generate_mutations(sequences, n=10, inversions=False, ksize=31, seed=None):
    if not seed:
        seed = random.randrange(sys.maxsize)
    print('[kevlar::gentrio] using random seed', seed, file=sys.stderr)
    rng = random.Random(seed)

    types = ['snv', 'ins', 'del']
    if inversions:
        # types.append('inv')
        raise NotImplementedError('feature pending')
    for _ in range(n):
        seqid = rng.choice(list(sequences.keys()))
        seq = sequences[seqid]
        seqlength = len(sequences[seqid])
        position = rng.randint(0, seqlength)
        muttype = rng.choice(types)

        if muttype == 'snv':
            offset = rng.randint(1, 3)
            refrseq, altseq, refrwindow, altwindow = mutate_snv(
                seq, position, offset, ksize
            )
        elif muttype == 'ins':
            length = rng.randint(5, 350)
            duplpos = rng.randint(0, seqlength)
            refrseq, altseq, refrwindow, altwindow = mutate_insertion(
                seq, position, length, duplpos, ksize
            )
        elif muttype == 'del':
            length = rng.randint(5, 350)
            refrseq, altseq, refrwindow, altwindow = mutate_deletion(
                seq, position, length, ksize
            )
        yield seqid, position, refrseq, altseq, refrwindow, altwindow
