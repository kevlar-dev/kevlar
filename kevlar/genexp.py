#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict, namedtuple
import random
import sys
import scipy.stats
import kevlar
from kevlar import MutableString
from kevlar.call import Variant


"""
- we need:
    - n common variants present in all samples (at some probability)
    - n' unique variants per variant (n' ~ Poisson(mu=2))
    - n'' case variants shared among all case samples (at some probability)
"""


def simulate_variants(sequences, numcommon=10, numcase=5, casesamples=10,
                      ctrlsamples=10, commonp=0.9, mu=2, casep=0.8,
                      weights=DWEIGHTS, rng=None):
    # Common variants
    mutator = generate_mutations(sequences, n=numcommon, weights=weights,
                                 rng=rng)
    for variant in mutator:
        casevector, ctrllvector = list(), list()
        while len(casevector) < casesamples:
            genotype = '1' if rng.random() < commonp else '0'
            casevector.append(genotype)
        while len(ctrllvector) < ctrlsamples:
            genotype = '1' if rng.random() < commonp else '0'
            ctrllvector.append(genotype)
        variant.info['CSV'] = ','.join(casevector)
        variant.info['CRV'] = ','.join(ctrllvector)
        yield variant

    # Unique variants
    for i in range(casesamples):
        casevector = ['0'] * casesamples
        ctrlvector = ['0'] * ctrlsamples
        casevector[i] = '1'
        numunique = scipy.stats.poisson.rvs(mu=2)
        mutator = generate_mutations(sequences, n=numunique, weights=weights,
                                     rng=rng)
        for variant in mutator:
            variant.info['CSV'] = ','.join(casevector)
            variant.info['CRV'] = ','.join(ctrlvector)
            yield variant
    for i in range(ctrlamples):
        casevector = ['0'] * casesamples
        ctrlvector = ['0'] * ctrlsamples
        ctrlvector[i] = '1'
        numunique = scipy.stats.poisson.rvs(mu=2)
        mutator = generate_mutations(sequences, n=numunique, weights=weights,
                                     rng=rng)
        for variant in mutator:
            variant.info['CSV'] = ','.join(casevector)
            variant.info['CRV'] = ','.join(ctrlvector)
            yield variant

    # Case variants
    mutator = generate_mutations(sequences, n=numcommon, weights=weights,
                                 rng=rng)
    ctrlvector = ['0'] * ctrlsamples
    for variant in mutator:
        casevector = list()
        while len(casevector) < casesamples:
            genotype = '1' if rng.random() < casep else '0'
            casevector.append(genotype)
        variant.info['CSV'] = ','.join(casevector)
        variant.info['CRV'] = ','.join(ctrllvector)
        yield variant


def genexp(sequences, outstreams, numcommon=20, numcase=10, weights=DWEIGHTS,
            seed=None, upint=100, logstream=sys.stderr):
    mutator = simulate_variants(
        sequences, numcommon=numcommon, numcase=numcase,
        casesamples=casesamples, ctrlsamples=ctrlsamples, commonp=commonp,
        mu=mu, casep=casep, weights=weights, rng=seed
    )
    variants = list(mutator)
    variants.sort(key=lambda v: v.position, reverse=True)

#    for seqid, sequence in sequences.items():
#        for ind in range(3):  # proband mother father
#            haploseqs = [MutableString(sequence), MutableString(sequence)]
#            for n, variant in enumerate(variants, 1):
#                if n % upint == 0:  # pragma: no cover
#                    message = (
#                        '    sequence={seq} individual={ind} '
#                        'variantsprocessed={vp}'.format(
#                            seq=seqid, ind=ind, vp=n
#                        )
#                    )
#                    print(message, file=logstream, flush=True)
#                if variant.seqid != seqid:
#                    continue
#                genotype = variant.genotypes[ind]
#                haplotypes = (genotype[0], genotype[2])
#                for hapindex in range(2):
#                    if haplotypes[hapindex] == '0':
#                        continue
#                    apply_mutation(
#                        haploseqs[hapindex], variant.position, variant._refr,
#                        variant._alt
#                    )
#            print('>', seqid, '_haplo1\n', haploseqs[0], sep='',
#                  file=outstreams[ind])
#            print('>', seqid, '_haplo2\n', haploseqs[1], sep='',
#                  file=outstreams[ind])
#
#    variants.sort(key=lambda v: (v.seqid, v.position))
#    for variant in variants:
#        yield variant


def main(args):
    pass
