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
from kevlar import MutableString
from kevlar.vcf import Variant


# Mappings for SNVs
nucl_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
index_to_nucl = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

# Default weights/probabilities for different mutation types
DWEIGHTS = {'snv': 0.8, 'ins': 0.1, 'del': 0.1}


"""
Viable inheritance scenarios

Each scenario is a 3-tuple showing the genotype of the trio for a particular
variant (child, mother, father).
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


def mutagenize(sequence, rng=None, rate=0.05):
    mutseq = list()
    for nucl in sequence:
        if rng and rng.random() < rate:
            offset = rng.choice([1, 2, 3])
            newindex = (nucl_to_index[nucl] + offset) % 4
            nucl = index_to_nucl[newindex]
        mutseq.append(nucl)
    return ''.join(mutseq)


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
        sequence[position+1:windowend]
    )

    return orignucl, newnucl, refrwindow, altwindow


def mutate_insertion(sequence, position, length, duplpos, rng=None, ksize=31):
    duplseq = mutagenize(sequence[duplpos:duplpos + length], rng, rate=0.05)
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
    altseq = sequence[position-1]
    refrseq = altseq + delseq

    windowstart = max(position - ksize + 1, 0)
    windowend = min(position + length + ksize - 1, len(sequence))
    refrwindow = sequence[windowstart:windowend]
    altwindow = '{:s}{:s}'.format(
        sequence[windowstart:position], sequence[position + length:windowend]
    )

    return refrseq, altseq, refrwindow, altwindow


def generate_mutations(sequences, n=10, ksize=31, weights=DWEIGHTS, rng=None):
    if rng is None:
        seed = random.randrange(sys.maxsize)
        print('[kevlar::gentrio] using random seed', seed, file=sys.stderr)
        rng = random.Random(seed)
    if isinstance(rng, int):
        rng = random.Random(rng)

    weightkeys = sorted(weights.keys())
    weightvalues = [weights[k] for k in weightkeys]
    for _ in range(n):
        seqid = rng.choice(list(sorted(sequences.keys())))
        seq = sequences[seqid]
        seqlength = len(sequences[seqid])
        position = rng.randint(0, seqlength - 1)
        muttype = weighted_choice(weightkeys, weightvalues, rng)

        if muttype == 'snv':
            offset = rng.randint(1, 3)
            refrseq, altseq, refrwindow, altwindow = mutate_snv(
                seq, position, offset, ksize
            )
        elif muttype == 'ins':
            length = rng.randint(5, 350)
            duplpos = rng.randint(0, seqlength)
            refrseq, altseq, refrwindow, altwindow = mutate_insertion(
                seq, position, length, duplpos, rng, ksize
            )
        elif muttype == 'del':
            length = rng.randint(5, 350)
            refrseq, altseq, refrwindow, altwindow = mutate_deletion(
                seq, position, length, ksize
            )
        else:
            raise ValueError('unknown mutation type {}'.format(muttype))
        yield Variant(seqid, position, refrseq, altseq, ALTWINDOW=altwindow,
                      REFRWINDOW=refrwindow)


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


def simulate_variant_genotypes(sequences, ninh=20, ndenovo=10,
                               weights=DWEIGHTS, rng=None):
    if rng is None:
        seed = random.randrange(sys.maxsize)
        print('[kevlar::gentrio] using random seed', seed, file=sys.stderr)
        rng = random.Random(seed)
    if isinstance(rng, int):
        rng = random.Random(rng)

    mutator = generate_mutations(sequences, n=ninh, weights=weights, rng=rng)
    for variant in mutator:
        genotypes = pick_inheritance_genotypes(rng)
        gtstring = ','.join(genotypes)
        variant.annotate('GT', gtstring)
        yield variant

    mut8r = generate_mutations(sequences, n=ndenovo, weights=weights, rng=rng)
    for variant in mut8r:
        genotypes = (rng.choice(['0/1', '1/0']), '0/0', '0/0')
        gtstring = ','.join(genotypes)
        variant.annotate('GT', gtstring)
        yield variant


def apply_mutation(sequence, position, refr, alt):
    if len(refr) == len(alt):  # SNV
        assert sequence[position] == refr
        sequence[position] = alt
    elif len(refr) < len(alt):  # Insertion
        sequence[position:position] = alt[1:]
    else:  # Deletion
        dellength = len(refr) - len(alt)
        del sequence[position:position+dellength]


def weights_str_to_dict(wstring):
    weights = dict()
    for keyvaluepair in wstring.split(','):
        muttype, relfreq = keyvaluepair.split('=')
        relfreq = float(relfreq)
        weights[muttype] = relfreq

    # Normalize!
    total = sum(weights.values())
    weights = {t: (v / total) for t, v in weights.items()}

    return weights


def gentrio(sequences, outstreams, ninh=20, ndenovo=10, weights=DWEIGHTS,
            seed=None, upint=100, logstream=sys.stderr):
    assert len(outstreams) == 3
    mutator = simulate_variant_genotypes(
        sequences, ninh=ninh, ndenovo=ndenovo, weights=weights, rng=seed
    )
    variants = list(mutator)
    variants.sort(key=lambda v: v.position, reverse=True)

    for seqid, sequence in sequences.items():
        for ind in range(3):  # proband mother father
            haploseqs = [MutableString(sequence), MutableString(sequence)]
            for n, variant in enumerate(variants, 1):
                if n % upint == 0:  # pragma: no cover
                    message = (
                        '    sequence={seq} individual={ind} '
                        'variantsprocessed={vp}'.format(
                            seq=seqid, ind=ind, vp=n
                        )
                    )
                    print(message, file=logstream, flush=True)
                if variant.seqid != seqid:
                    continue
                genotype = variant.genotypes[ind]
                haplotypes = (genotype[0], genotype[2])
                for hapindex in range(2):
                    if haplotypes[hapindex] == '0':
                        continue
                    apply_mutation(
                        haploseqs[hapindex], variant.position, variant._refr,
                        variant._alt
                    )
            print('>', seqid, '_haplo1\n', haploseqs[0], sep='',
                  file=outstreams[ind])
            print('>', seqid, '_haplo2\n', haploseqs[1], sep='',
                  file=outstreams[ind])

    variants.sort(key=lambda v: (v.seqid, v.position))
    for variant in variants:
        yield variant


def main(args):
    timer = kevlar.Timer()
    timer.start()

    timer.start('loadgenome')
    print('[kevlar::gentrio] Loading genome...', end='', file=sys.stderr)
    seqfile = kevlar.open(args.genome, 'r')
    genomeseqs = kevlar.seqio.parse_seq_dict(seqfile)
    elapsed = timer.stop('loadgenome')
    print('done! ({:.3f} seconds elapsed)'.format(elapsed), file=sys.stderr)

    samples = ('proband', 'mother', 'father')
    outfiles = ['{:s}-{:s}.fasta'.format(args.prefix, s) for s in samples]
    outstreams = [kevlar.open(outfile, 'w') for outfile in outfiles]

    vcfout = None
    if args.vcf:
        vcfout = kevlar.open(args.vcf, 'w')
        kevlar.vcf_header(vcfout, source='kevlar::gentrio', infoheader=True)
    weights = weights_str_to_dict(args.weights)
    mutator = gentrio(
        genomeseqs, outstreams, ninh=args.inherited, ndenovo=args.de_novo,
        weights=weights, seed=args.seed, logstream=args.logfile
    )

    timer.start('mutate')
    print('[kevlar::gentrio] Begin generating and applying mutations:',
          file=sys.stderr)
    for variant in mutator:
        if vcfout:
            print(variant.vcf, file=vcfout)
    elapsed = timer.stop('mutate')
    print('[kevlar::gentrio] Done applying mutations! ', end='',
          file=sys.stderr)
    print('({:.3f} seconds elapsed)'.format(elapsed), file=sys.stderr)

    for outstream in outstreams:
        outstream.close()

    elapsed = timer.stop()
    print('[kevlar::gentrio] Trio simulation complete; ', file=sys.stderr)
    print(' total runtime: {:.3f} seconds'.format(elapsed), file=sys.stderr)
