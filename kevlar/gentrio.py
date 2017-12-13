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

"""
Viable inheritance scenarios

Each scenario is a 3-tuple with showing the genotype of the trio for a
particular variant (child, mother, father).
- 0 refers to homozygous reference (0/0)
- 1 refers to heterzygous (0/1 or 1/0)
- 2 refers to homozygous alternate (1/1)

Combinations of genotypes that are not valid inheritance scenarios are
excluded.
"""
inheritance_scenarios = [

    (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 1), (1, 0, 2), (1, 1, 0),
    (1, 1, 1), (1, 1, 2), (1, 2, 0), (1, 2, 1), (2, 1, 1), (2, 1, 2),
    (2, 2, 1), (2, 2, 2),
]


def weighted_choice(values, weights, rng=random.Random()):
    """Stolen shamelessly from https://stackoverflow.com/a/3679747/459780."""
    assert len(values) == len(weights)
    total = sum(weights)
    r = rng.uniform(0, total)
    cumsum = 0
    for v, w in zip(values, weights):
        if cumsum + w >= r:
            return v
        cumsum += w
    assert False  # pragma: no cover


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


def generate_mutations(sequences, n=10, inversions=False, ksize=31, rng=None):
    if not rng:
        seed = random.randrange(sys.maxsize)
        print('[kevlar::gentrio] using random seed', seed, file=sys.stderr)
        rng = random.Random(seed)
    elif isinstance(rng, int):
        rng = random.Random(rng)

    types = ['snv', 'ins', 'del']
    weights = [0.7, 0.15, 0.15]
    if inversions:
        # types.append('inv')
        raise NotImplementedError('feature pending')

    for _ in range(n):
        seqid = rng.choice(list(sequences.keys()))
        seq = sequences[seqid]
        seqlength = len(sequences[seqid])
        position = rng.randint(0, seqlength)
        muttype = weighted_choice(types, weights, rng)

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
        else:
            assert muttype == 'del'
            length = rng.randint(5, 350)
            refrseq, altseq, refrwindow, altwindow = mutate_deletion(
                seq, position, length, ksize
            )
        yield seqid, position, refrseq, altseq, refrwindow, altwindow


def pick_inheritance_genotypes(rng):
    genotype_codes = rng.choice(inheritance_scenarios)
    genotypes = list()
    for code in genotype_codes:
        if code == 0:
            genotype = '0/0'
        elif code == 2:
            genotype = '1/1'
        else:
            genotype = rng.choice(['0/1', '1/0'])
        genotypes.append(genotype)
    return tuple(genotypes)


def simulate_variant_genotypes(sequences, ninh=20, ndenovo=10, seed=None):
    if not seed:
        seed = random.randrange(sys.maxsize)
    print('[kevlar::gentrio] using random seed', seed, file=sys.stderr)
    rng = random.Random(seed)

    variants = list()
    mutator = generate_mutations(sequences, n=ninh, rng=rng)
    for variant in mutator:
        genotypes = pick_inheritance_genotypes(rng)
        variants.append((variant, genotypes))

    mutator = generate_mutations(sequences, n=ndenovo, rng=rng)
    for variant in mutator:
        genotypes = (rng.choice(['0/1', '1/0']), '0/0', '0/0')
        variants.append((variant, genotypes))

    variants.sort(key=lambda v: (v[0][0], v[0][1]))
    for variant, genotypes in variants:
        yield variant, genotypes
