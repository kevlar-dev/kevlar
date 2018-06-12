#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
import json
import math
import sys
import threading

import kevlar
import khmer
import numpy
import pandas


class KevlarZeroAbundanceDistError(ValueError):
    pass


def count_first_pass(infiles, counts, mask, nthreads=1, logstream=sys.stderr):
    message = 'Processing input with {:d} threads'.format(nthreads)
    print('[kevlar::dist]', message, file=logstream)

    for filename in infiles:
        print('    -', filename, file=logstream)
        parser = khmer.ReadParser(filename)
        threads = list()
        for _ in range(nthreads):
            thread = threading.Thread(
                target=counts.consume_seqfile_with_mask,
                args=(parser, mask,),
                kwargs={'threshold': 1, 'consume_masked': True},
            )
            threads.append(thread)
            thread.start()
        for thread in threads:
            thread.join()

    print('[kevlar::dist] Done processing input!', file=logstream)


def count_second_pass(infiles, counts, nthreads=1, logstream=sys.stderr):
    print('[kevlar::dist] Second pass over the data', file=logstream)
    tracking = khmer.Nodetable(counts.ksize(), 1, 1, primes=counts.hashsizes())
    abund_lists = list()

    def __do_abund_dist(parser):
        abund = counts.abundance_distribution(parser, tracking)
        abund_lists.append(abund)

    for filename in infiles:
        print('    -', filename, file=logstream)
        parser = khmer.ReadParser(filename)
        threads = list()
        for _ in range(nthreads):
            thread = threading.Thread(
                target=__do_abund_dist,
                args=(parser,)
            )
            threads.append(thread)
            thread.start()
        for thread in threads:
            thread.join()

    assert len(abund_lists) == len(infiles) * nthreads
    abundance = defaultdict(int)
    for abund in abund_lists:
        for i, count in enumerate(abund):
            if i > 0 and count > 0:
                abundance[i] += count

    print('[kevlar::dist] Done second pass over input!', file=logstream)

    return abundance


def weighted_mean_std_dev(values, weights):
    mu = numpy.average(values, weights=weights)
    sigma = math.sqrt(numpy.average((values-mu)**2, weights=weights))
    return mu, sigma


def calc_mu_sigma(abundance):
    total = sum(abundance.values())
    if total == 0:
        message = 'all k-mer abundances are 0, please check input files'
        raise KevlarZeroAbundanceDistError(message)
    mu, sigma = weighted_mean_std_dev(
        list(abundance.keys()),
        list(abundance.values()),
    )
    return mu, sigma


def compute_dist(abundance):
    total = sum(abundance.values())
    fields = ['Abundance', 'Count', 'CumulativeCount', 'CumulativeFraction']
    data = pandas.DataFrame(columns=fields)
    cuml = 0
    for abund, count in sorted(abundance.items()):
        assert count > 0, (abund, count)
        cuml += count
        frac = cuml / total
        row = {
            'Abundance': abund,
            'Count': count,
            'CumulativeCount': cuml,
            'CumulativeFraction': frac,
        }
        data = data.append(row, ignore_index=True)
    return data


def dist(infiles, mask, ksize=31, memory=1e6, threads=1, logstream=sys.stderr):
    counts = khmer.Counttable(ksize, memory / 4, 4)
    count_first_pass(infiles, counts, mask, nthreads=threads,
                     logstream=logstream)
    abundance = count_second_pass(infiles, counts, nthreads=threads,
                                  logstream=logstream)
    mu, sigma = calc_mu_sigma(abundance)
    data = compute_dist(abundance)
    return mu, sigma, data


def main(args):
    mask = khmer.Nodetable.load(args.mask)
    mu, sigma, data = dist(
        args.infiles, mask, ksize=args.ksize, memory=args.memory,
        threads=args.threads, logstream=args.logfile
    )
    out = {'mu': mu, 'sigma': sigma}
    print(json.dumps(out))

    if args.tsv:
        data.to_csv(args.tsv, sep='\t', index=False)

    if args.plot:
        import os
        try:  # pragma: no cover
            import matplotlib
            if os.environ.get('DISPLAY', '') == '':
                matplotlib.use('Agg')
            from matplotlib import pyplot as plt
        except RuntimeError as rerr:  # pragma: no cover
            if 'Python is not installed as a framework' not in str(rerr):
                raise rerr
            message = 'There was a problem loading matplotlib. '
            message += 'Try https://stackoverflow.com/q/21784641/459780 '
            message += 'for troubleshooting ideas.'
            raise RuntimeError(message)
        matplotlib.rcParams["figure.figsize"] = [12, 6]
        matplotlib.rcParams['axes.labelsize'] = 16
        matplotlib.rcParams['xtick.labelsize'] = 16
        plt.plot(data['Abundance'], data['Count'], color='blue')
        plt.axvline(x=mu, color='blue', linestyle='--')
        plt.axvline(x=mu - sigma, color='red', linestyle=':')
        plt.axvline(x=mu + sigma, color='red', linestyle=':')
        plt.xlim(args.plot_xlim)
        plt.xlabel('K-mer abundance')
        plt.ylabel('Frequency')
        plt.savefig(args.plot, dpi=300)
