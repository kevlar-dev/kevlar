#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
import math
import sys
import threading

import kevlar
import khmer
import numpy


class KevlarZeroAbundanceDistError(ValueError):
    pass


def weighted_mean_std_dev(values, weights):
    mu = numpy.average(values, weights=weights)
    sigma = math.sqrt(numpy.average((values-mu)**2, weights=weights))
    return mu, sigma


def count_first_pass(infiles, counts, mask, threads=1, logstream=sys.stderr):
    message = 'Processing input with {:d} threads'.format(threads)
    print('[kevlar::dist]', message, file=logstream)

    for filename in infiles:
        print('    -', filename, file=logstream)
        parser = khmer.ReadParser(filename)
        threads = list()
        for _ in range(threads):
            thread = threading.Thread(
                target=counts.consume_seqfile_with_mask,
                args=(parser, mask,),
                kwargs={'threshold': 1, 'complement': True},
            )
            threads.append(thread)
            thread.start()
        for thread in threads:
            thread.join()

    print('[kevlar::dist] Done processing input!', file=logstream)


def count_second_pass(infiles, counts, threads=1, logstream=logstream):
    print('[kevlar::dist] Second pass over the data', file=logstream)
    tracking = khmer.Nodetable(ksize, 1, 1, primes=counttable.hashsizes())
    abund_lists = list()

    def __do_abund_dist(parser):
        abund = counts.abundance_distribution(parser, tracking)
        abund_lists.append(abund)

    for filename in infiles:
        print('    -', filename, file=logstream)
        parser = khmer.ReadParser(filename)
        threads = list()
        for _ in range(threads):
            thread = threading.Thread(
                target=__do_abund_dist,
                args=(parser,)
            )
            threads.append(thread)
            thread.start()
        for thread in threads:
            thread.join()

    assert len(abund_lists) == len(infiles) * threads
    abundance = defaultdict(int)
    for abund in abund_lists:
        for i, count in enumerate(abund):
            abundance[i] += count

    print('[kevlar::dist] Done second pass over input!', file=logstream)

    return abundance


def calc_mu_sigma(abundance):
    message = 'ignoring {:d} k-mers as not in mask'.format(abundance[0])
    print('[kevlar::dist]', message, file=logstream)
    del(abundance[0])
    total = sum(abundance.values())
    if total == 0:
        message = 'all k-mer abundances are 0, please check input files'
        raise KevlarZeroAbundanceDistError(message)
    mu, sigma = weighted_mean_std_dev(abundance.keys(), abundance.values())
    return mu_sigma


def compute_dist(abundance):
    total = sum(abundance.values())
    fields = ['Abundance', 'Count', 'CumulativeCount', 'CumulativeFraction']
    data = pandas.DataFrame(columns=fields)
    cuml = 0
    for abund, count in sorted(abundance.items()):
        if count == 0:
            continue
        cuml += count
        frac = cuml / total
        row = {
            'Abundance': abund,
            'Count': count,
            'CumulativeCount': cuml,
            'CumulativeFraction': frac,
        }
        data = data.append(row, ignore_index=True)
    return abund


def dist(infiles, mask, ksize=31, memory=1e6, threads=1, logstream=sys.stderr):
    counts = khmer.Counttable(ksize, memory / 4, 4)
    count_first_pass(infiles, counts, mask, threads=threads,
                     logstream=logstream)
    abundance = count_second_pass(infiles, counts, threads=threads,
                                  logstream=logstream)
    mu, sigma = calc_mu_sigma(abundance)
    data = compute_dist(abundance)
    return mu, sigma, data
