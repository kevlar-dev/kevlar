#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import kevlar
from kevlar.vcf import Variant
from math import log
import scipy.stats


def get_abundances(altseq, refrseq, case, controls, refr):
    """Create a nested list of k-mer abundances.

    abundances = [
        [15, 14, 13, 16, 14, 15, 14, 14],  # k-mer abundances from case/proband
        [0, 0, 1, 0, 2, 10, 0, 0],  # k-mer abundances from parent/control 1
        [0, 1, 1, 0, 1, 0, 2, 0],  # k-mer abundances from parent/control 2
    ]
    refr_abunds = [1, 1, 2, 1, 4, 2, 1, 1]  # genomic freq of refr allele kmers
    """
    altkmers = case.get_kmers(altseq)
    refrkmers = case.get_kmers(refrseq)
    if len(altseq) == len(refrseq):
        valid_alt_kmers = list()
        refr_abunds = list()
        for altkmer, refrkmer in zip(altkmers, refrkmers):
            if refr.get(altkmer) == 0:
                valid_alt_kmers.append(altkmer)
                refr_abunds.append(refr.get(refrkmer))
    else:
        valid_alt_kmers = [k for k in altkmers if refr.get(k) == 0]
        refr_abunds = [None] * len(valid_alt_kmers)
    ndropped = len(altkmers) - len(valid_alt_kmers)

    abundances = list()
    for _ in range(len(controls) + 1):
        abundances.append(list())
    for kmer in valid_alt_kmers:
        abund = case.get(kmer)
        abundances[0].append(abund)
        for control, abundlist in zip(controls, abundances[1:]):
            abund = control.get(kmer)
            abundlist.append(abund)
    return abundances, refr_abunds, ndropped


def abund_log_prob(genotype, abundance, refrabund=None, mean=30.0, sd=8.0,
                   error=0.001):
    """Calculate probability of k-mer abundance conditioned on genotype.

    The `genotype` variable represents the number of assumed allele copies and
    is one of {0, 1, 2} (corresponding to genotypes {0/0, 0/1, and 1/1}). The
    `mean` and `sd` variables describe a normal distribution of observed
    abundances of k-mers with copy number 2. The `error` parameter is the error
    rate.

    For SNVs, there is a 1-to-1 correspondence of alternate allele k-mers to
    reference allele k-mers. We can therefore check the frequency of the
    reference allele in the reference genome and increase the error rate if it
    is repetitive. There is no such mapping of alt allele k-mers to refr allele
    k-mers for indels, so we use a flat error rate.
    """
    if genotype == 0:
        if refrabund:  # SNV
            erate = error * mean * refrabund
        else:  # INDEL
            erate = error * mean * 0.1
        return abundance * log(erate)
    elif genotype == 1:
        return scipy.stats.norm.logpdf(abundance, mean / 2, sd / 2)
    elif genotype == 2:
        return scipy.stats.norm.logpdf(abundance, mean, sd)


def likelihood_denovo(abunds, refrabunds, mean=30.0, sd=8.0, error=0.001):
    assert len(abunds[1]) == len(refrabunds)
    assert len(abunds[2]) == len(refrabunds)
    logsum = 0.0

    # Case
    for abund in abunds[0]:
        logsum += abund_log_prob(1, abund, mean=mean, sd=sd)
    # Controls
    for altabunds in abunds[1:]:
        for alt, refr in zip(altabunds, refrabunds):
            logsum += abund_log_prob(0, alt, refrabund=refr, mean=mean,
                                     error=error)
    return logsum


def likelihood_false(abunds, refrabunds, mean=30.0, error=0.001):
    assert len(abunds[1]) == len(refrabunds)
    assert len(abunds[2]) == len(refrabunds)
    logsum = 0.0
    for altabunds in abunds:
        for alt, refr in zip(altabunds, refrabunds):
            logsum += abund_log_prob(0, alt, refrabund=refr, mean=mean,
                                     error=error)
    return logsum


def likelihood_inherited(abunds, mean=30.0, sd=8.0, error=0.001):
    """Compute the likelihood that a variant is inherited.

    There are 15 valid inheritance scenarios, 11 of which (shown below) result
    in the proband carrying at least one copy of the alternate allele. Select
    the one with the highest likelihood.

    The other likelihood calculations are implemented to handle an arbitrary
    number of controls, but this can only handle trios.
    """
    scenarios = [
        (1, 0, 1), (1, 0, 2),
        (1, 1, 0), (1, 1, 1), (1, 1, 2),
        (1, 2, 0), (1, 2, 1),
        (2, 1, 1), (2, 1, 2),
        (2, 2, 1), (2, 2, 2),
    ]
    logsum = 0.0
    abundances = zip(abunds[0], abunds[1], abunds[2])
    for a_c, a_m, a_f in abundances:
        maxval = None
        for g_c, g_m, g_f in scenarios:
            p_c = abund_log_prob(g_c, a_c, mean=mean, sd=sd, error=error)
            p_m = abund_log_prob(g_m, a_m, mean=mean, sd=sd, error=error)
            p_f = abund_log_prob(g_f, a_f, mean=mean, sd=sd, error=error)
            testsum = p_c + p_m + p_f + log(1.0 / 15.0)
            if maxval is None or testsum > maxval:
                maxval = testsum
        logsum += maxval
    return log(15.0 / 11.0) + logsum  # 1 / (11/15)
