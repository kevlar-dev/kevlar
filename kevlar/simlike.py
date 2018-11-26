#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
from math import log
import sys

import kevlar
from kevlar.vcf import Variant
import khmer
import scipy.stats


class KevlarSampleLabelingError(ValueError):
    pass


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
                   error=0.001, dynamic=True):
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
            erate = error * mean
            if dynamic:
                erate *= refrabund
        else:  # INDEL
            erate = error * mean * 0.1
        return abundance * log(erate)
    elif genotype == 1:
        return scipy.stats.norm.logpdf(abundance, mean / 2, sd / 2)
    elif genotype == 2:
        return scipy.stats.norm.logpdf(abundance, mean, sd)


def likelihood_denovo(abunds, refrabunds, mean=30.0, sd=8.0, error=0.001,
                      dynamic=True):
    assert len(abunds[1]) == len(refrabunds)
    assert len(abunds[2]) == len(refrabunds)
    logsum = 0.0

    # Case
    for abund in abunds[0]:
        logsum += abund_log_prob(1, abund, mean=mean, sd=sd, dynamic=dynamic)
    # Controls
    for altabunds in abunds[1:]:
        for alt, refr in zip(altabunds, refrabunds):
            logsum += abund_log_prob(0, alt, refrabund=refr, mean=mean,
                                     error=error, dynamic=dynamic)
    return logsum


def likelihood_false(abunds, refrabunds, mean=30.0, error=0.001,
                     dynamic=True):
    assert len(abunds[1]) == len(refrabunds)
    assert len(abunds[2]) == len(refrabunds)
    logsum = 0.0
    for altabunds in abunds:
        for alt, refr in zip(altabunds, refrabunds):
            logsum += abund_log_prob(0, alt, refrabund=refr, mean=mean,
                                     error=error, dynamic=dynamic)
    return logsum


def likelihood_inherited(abunds, mean=30.0, sd=8.0, error=0.001,
                         dynamic=True):
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
            p_c = abund_log_prob(g_c, a_c, mean=mean, sd=sd, error=error,
                                 dynamic=dynamic)
            p_m = abund_log_prob(g_m, a_m, mean=mean, sd=sd, error=error,
                                 dynamic=dynamic)
            p_f = abund_log_prob(g_f, a_f, mean=mean, sd=sd, error=error,
                                 dynamic=dynamic)
            testsum = p_c + p_m + p_f + log(1.0 / 15.0)
            if maxval is None or testsum > maxval:
                maxval = testsum
        logsum += maxval
    return log(15.0 / 11.0) + logsum  # 1 / (11/15)


def joinlist(thelist):
    if len(thelist) == 0:
        return '.'
    else:
        return ','.join([str(v) for v in thelist])


def calc_likescore(call, altabund, refrabund, mu, sigma, epsilon,
                   dynamic=True):
    lldn = likelihood_denovo(altabund, refrabund, mean=mu, sd=sigma,
                             error=epsilon, dynamic=dynamic)
    llfp = likelihood_false(altabund, refrabund, mean=mu, error=epsilon,
                            dynamic=dynamic)
    llih = likelihood_inherited(altabund, mean=mu, sd=sigma, error=epsilon,
                                dynamic=dynamic)
    likescore = lldn - max(llfp, llih)
    call.annotate('LLDN', lldn)
    call.annotate('LLFP', llfp)
    call.annotate('LLIH', llih)
    call.annotate('LIKESCORE', likescore)


def default_sample_labels(nsamples):
    samples = list()
    for i in range(nsamples):
        samples.append('Control{:d}'.format(i))
    samples[0] = 'Case'
    return samples


def annotate_abundances(call, abundances, samplelabels):
    for sample, abundlist in zip(samplelabels, abundances):
        abundstr = joinlist(abundlist)
        call.format(sample, 'ALTABUND', abundstr)


def process_partition(partitionid, calls):
    maxscore = max([c.attribute('LIKESCORE') for c in calls])
    maxcalls = list()
    for call in calls:
        if call.attribute('LIKESCORE') == maxscore:
            maxcalls.append(call)
        else:
            call.filter(kevlar.vcf.VariantFilter.PartitionScore)
    for call in maxcalls:
        call.annotate('CALLCLASS', partitionid)
    if calls[0].attribute('MATEDIST') is not None:
        matedists = set([c.attribute('MATEDIST') for c in calls])
        if matedists == set([float('inf')]):
            for call in calls:
                call.filter(kevlar.vcf.VariantFilter.MateFail)


def simlike(variants, case, controls, refr, mu=30.0, sigma=8.0, epsilon=0.001,
            dynamic=True, casemin=5, samplelabels=None, logstream=sys.stderr):
    calls_by_partition = defaultdict(list)
    if samplelabels is None:
        samplelabels = default_sample_labels(len(controls) + 1)
    for call in variants:
        if call.window is None or len(call.window) < case.ksize():
            if call.filterstr == 'PASS':
                msg = 'WARNING: stubbornly refusing to compute likelihood for '
                if call.window is None:
                    msg += 'variant with no spanning window'
                else:
                    msg += 'variant-spanning window {:s}'.format(call.window)
                    msg += ', shorter than k size {:d}'.format(case.ksize())
                print('[kevlar::simlike]', msg, file=logstream)
            call.annotate('LIKESCORE', float('-inf'))
            calls_by_partition[call.attribute('PART')].append(call)
            continue
        altabund, refrabund, ndropped = get_abundances(
            call.window, call.refrwindow, case, controls, refr
        )
        call.annotate('DROPPED', ndropped)
        abovethresh = [a for a in altabund[0] if a > casemin]
        if len(abovethresh) == 0:
            call.filter(kevlar.vcf.VariantFilter.PassengerVariant)
        calc_likescore(call, altabund, refrabund, mu, sigma, epsilon,
                       dynamic=dynamic)
        annotate_abundances(call, altabund, samplelabels)
        calls_by_partition[call.attribute('PART')].append(call)

    allcalls = list()
    for partition, calls in calls_by_partition.items():
        process_partition(partition, calls)
        allcalls.extend(calls)

    allcalls.sort(key=lambda c: c.attribute('LIKESCORE'), reverse=True)
    for call in allcalls:
        if call.attribute('LIKESCORE') < 0.0:
            call.filter(kevlar.vcf.VariantFilter.LikelihoodFail)
        yield call


def main(args):
    nsamples = len(args.controls) + 1
    if args.sample_labels:
        nlabels = len(args.sample_labels)
        if nlabels and nlabels != nsamples:
            message = 'provided {:d} labels'.format(nlabels)
            message += ' but {:d} samples'.format(nsamples)
            raise KevlarSampleLabelingError(message)
    else:
        args.sample_labels = default_sample_labels(nsamples)

    case = khmer.Counttable.load(args.case)
    controls = [khmer.Counttable.load(c) for c in args.controls]
    refr = khmer.SmallCounttable.load(args.refr)

    instream = kevlar.open(args.vcf, 'r')
    outstream = kevlar.open(args.out, 'w')
    reader = kevlar.vcf.VCFReader(instream)
    writer = kevlar.vcf.VCFWriter(outstream, source='kevlar::simlike')

    for label in args.sample_labels:
        writer.register_sample(label)
    writer.write_header()

    calculator = simlike(
        reader, case, controls, refr, mu=args.mu, sigma=args.sigma,
        epsilon=args.epsilon, dynamic=args.dynamic, casemin=args.case_min,
        samplelabels=args.sample_labels, logstream=args.logfile,
    )
    for call in calculator:
        writer.write(call)
