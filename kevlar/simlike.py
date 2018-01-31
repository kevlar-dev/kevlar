#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import khmer
import kevlar
from kevlar.vcf import Variant
import scipy


def KevlarUnknownGenotypeError(ValueError):
    pass


def filter_refr(akmers, rkmers, refrcounts):
    """
    Remove k-mers that are present in the reference.

    The `akmers` and `rkmers` variables are lists of k-mers that contain an
    alternate and reference allele, respectively. The `refr` variable is a
    table of k-mer counts from the reference genome. Any alternate allele k-mer
    is discarded if it is present in the reference genome, and any reference
    allele k-mer is discarded if it is present in multiple copies in the
    reference genome.
    """
    newa, newr = list(), list()
    for a in akmers:
        if refr.get(a) == 0:
            newa.append(a)
    for r in rkmers:
        if refr.get(r) < 2:
            newr.append(r)
    return newa, newr


def get_abundances(kmers, casecounts, controlcounts):
    """
    Create a nested list of k-mer abundances.

    Given a list of k-mers, store the k-mer abundances using one sub-list per
    sample, as shown below.

    abunds = [
        [15, 14, 13, 16, 14, 15, 14, 14],  # k-mer abundances from case/proband
        [0, 0, 1, 0, 2, 10, 0, 0],  # k-mer abundances from parent/control 1
        [0, 1, 1, 0, 1, 0, 2, 0],  # k-mer abundances from parent/control 2
    ]
    """
    nsamples = len(controlcounts) + 1
    abundances = [list()] * nsamples
    for kmer in kmers:
        abund = casecounts.get(kmer)
        abundances[0].append(abund)
        for i in range(len(controlcounts)):
            abund = controlcounts[i].get(kmer)
            abundances[i+1].append(abund)
    return abundances


def set_error_rates(error, nsamples):
    """Set error rate dynamically"""
    if isinstance(error, float):
        errors = [error] * nsamples
    elif isinstance(error, list):
        assert len(error) == nsamples
        for e in error:
            assert isinstance(e, float)
        errors = error
    else:
        message = 'variable {}'.format(error)
        message += 'doesn\'t look/smell/quack like a float or list of floats'
        raise ValueError(message)
    return errors


def abund_log_prob(genotype, abundance, mean=40.0, sd=8.0, error=0.01):
    """
    Calculate conditional k-mer abundance probability

    Compute the log probability of the given k-mer abundance conditioned on the
    given genotype (copy number). The `genotype` variable represents the number
    of assumed allele copies and is one of {0, 1, 2} (corresponding to
    genotypes {0/0, 0/1, and 1/1}). The `mean` and `sd` variables describe a
    normal distribution of observed abundances of k-mers with copy number 2.
    The `error` parameter is the error rate.
    """
    if genotype not in (0, 1, 2):
        raise KevlarUnknownGenotypeError(gentype)
    if genotype == 0:
        return abundance * log(error)
    if genotype == 1:
        p = scipy.stats.norm.logpdf(abundance, mean / 2, sd / 2)
        return log(p)
    if genotype == 2:
        p = scipy.stats.norm.logpdf(abundance, mean, sd)
        return log(p)


def likelihood_denovo(altabunds, refrabunds, mean=30.0, sd=8.0, error=0.01):
    """
    Compute the likelihood that a variant is de novo.

    The `altabunds` and `refrabunds` variables should store nested list of
    k-mer abundances, structured as shown below.

    abunds = [
        [15, 14, 13, 16, 14, 15, 14, 14],  # k-mer abundances from proband
        [0, 0, 1, 0, 2, 10, 0, 0],  # k-mer abundances from parent/control 1
        [0, 1, 1, 0, 1, 0, 2, 0],  # k-mer abundances from parent/control 2
    ]

    The `error` variable can be a single float or a list of floats. If it is a
    single float, the same error rate will be applied to all control samples.
    If it is a list of floats, it must contain 1 value per control sample.
    """
    errors = set_error_rates(error, nsamples=len(altabunds))
    logsum = 0.0
    # Case
    for abund in altabunds[0] + refrabunds[0]:
        logsum += abund_log_prob(1, abund, mean=mean, sd=sd)
    # Controls
    for alt, refr, err in zip(altabunds[1:], refrabunds[1:], errors[1:]):
        for aabund in alt:
            logsum += abund_log_prob(0, aabund, error=err)
        for rabund in refr:
            logsum += abund_log_prob(2, rabund, mean=mean, sd=sd)
    return logsum


def likelihood_false(altabunds, error=0.01):
    """
    Compute the likelihood that a variant is false.

    The `altabunds` variable should store nested list of k-mer abundances,
    structured as shown below.

    abunds = [
        [15, 14, 13, 16, 14, 15, 14, 14],  # k-mer abundances from proband
        [0, 0, 1, 0, 2, 10, 0, 0],  # k-mer abundances from parent/control 1
        [0, 1, 1, 0, 1, 0, 2, 0],  # k-mer abundances from parent/control 2
    ]

    The `error` variable can be a single float or a list of floats. If it is a
    single float, the same error rate will be applied to all case/proband and
    control/parental samples. If it is a list of floats, it must contain 1
    value per sample.
    """
    errors = set_error_rates(error, nsamples=len(altabunds))
    logsum = 0.0
    for abundlist, e in zip(altabunds, errors):
        for abund in abundlist:
            logsum += abund_log_prob(0, abund, error=e)
    return logsum


def likelihood_inherited(altabunds, mean=30.0, sd=8.0, error=0.01):
    """
    Compute the likelihood that a variant is inherited.

    There are 15 valid inheritance scenarios, 11 of which (shown below) result
    in the proband carrying at least one copy of the alternate allele. Sum
    probabilities across all these different scenarios to derive an aggregate
    likelihood that the variant is inherited.

    The other likelihood calculations are implemented to handle an arbitrary
    number of controls, but this can only handle trios.
    """
    scenarios = (
        (1, 0, 1), (1, 0, 2),
        (1, 1, 0), (1, 1, 1), (1, 1, 2),
        (1, 2, 0), (1, 2, 1),
        (2, 1, 1), (2, 1, 2),
        (2, 2, 1), (2, 2, 2),
    )
    errors = set_error_rates(error, nsamples=3)
    logsum = 0.0
    abundances = zip(altabunds[0], altabunds[1], altabunds[2])
    for a_c, a_m, a_f in abundances:
        maxval = None
        for g_c, g_m, g_f in scenarios:
            testsum = abund_log_prob(g_c, a_c, mean, sd, errors[0]) + \
                      abund_log_prob(g_m, a_m, mean, sd, errors[1]) + \
                      abund_log_prob(g_f, a_f, mean, sd, errors[2]) + \
                      log(1.0 / 15.0)
            if maxval is None or testsum > maxval:
                maxval = testsum
        logsum += maxval
    return log(15.0 / 11.0) * logsum  # 1 / (11/15)


def compute_likelihoods(variant, casecounts, controlcounts, refr=None,
                        mean=30.0, sd=8.0, error=0.01, caselabel=None,
                        ctrllabels=None):
    if variant.window is None or variant.refrwindow is None:
        variant.info['DN'] = '-inf'
        return

    altkmers = casecounts.get_kmers(variant.window)
    refrkmers = casecounts.get_kmers(variant.refrwindow)
    altdrop, refrdrop = 0, 0
    if refr:
        nalt = len(altkmers)
        nrefr = len(refrkmers)
        altkmers, refrkmers = filter_refr(altkmers, refrkmers, refr)
        altdrop = nalt - len(altkmers)
        refrdrop = nrefr - len(refrkmers)
        if len(altkmers) == 0:
            variant.info['DN'] = '-inf'
            variant.info['NC'] = 'allaltinreference'
            # FIXME add new object attribute for filter pass or no
            return

    altabunds = get_abundances(altkmers, casecounts, controlcounts)
    refrabunds = get_abundances(refrkmers, casecounts, controlcounts)
    samplelabels = [caselabel] + ctrllabels
    for label, abundlist in zip(samplelabels, altabunds):
        abundstr = ','.join([str(a) for a in abundlist])
        variant.format(label, 'AA', abundstr)

    variant.info['FP'] = likelihood_false(altabunds, error=error)
    variant.info['IH'] = likelihood_inherited(altabunds, mean=mean, sd=sd,
                                              error=error)
    variant.info['DN'] = likelihood_denovo(altabunds, refrabunds, mean=mean,
                                           sd=sd, error=error)
    return
