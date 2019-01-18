#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019 Battelle National Biodefense Institute
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
import kevlar
from kevlar.intervalforest import IntervalForest


def populate_index_from_bed(instream):
    index = IntervalForest()
    for line in instream:
        if line.startswith('#') or line.strip() == '':
            continue
        values = line.strip().split()
        chrom = values[0]
        start, end = [int(coord) for coord in values[1:3]]
        strrepr = '{:s}:{:d}-{:d}'.format(chrom, start, end)
        index.insert(chrom, start, end, strrepr)
    return index


def compact(variants, index, delta=10):
    """Compact variants by call class

    Variant calls labeled with the same `CALLCLASS` attribute were predicted
    from the same partition (set of reads). While more than one of these may
    indeed be true variants with respect to the reference, we expect only one
    *de novo* variant per partition.

    This function assumes the variant calls are sorted by likelihood score (the
    `LIKESCORE` attribute in kevlar). Any call that does not pass filters is
    ignored. Then, for each CALLCLASS with multiple passing calls, all calls
    are discarded except for the one matching the true variant. If the
    CALLCLASS has no calls matching a true variant, all of the calls are
    discarded except for the highest scoring call.
    """
    variants_by_class = defaultdict(list)
    calls = list()
    for varcall in variants:
        if varcall.filterstr != 'PASS':
            continue
        callclass = varcall.attribute('CALLCLASS')
        if callclass is None:
            calls.append(varcall)
        else:
            variants_by_class[callclass].append(varcall)

    for callclass, calllist in variants_by_class.items():
        nmatches = 0
        match = None
        for varcall in calllist:
            hits = index.query(varcall.seqid, varcall.position, delta=delta)
            if hits == set():
                continue
            else:
                nmatches += 1
                if match is None:
                    match = varcall
        if nmatches == 0:
            calllist[0].annotate('EVAL', 'False')
            calls.append(calllist[0])
        else:
            assert nmatches > 0, nmatches
            if nmatches > 1:
                print('WARNING: found', nmatches, 'matches for CALLCLASS',
                      callclass, file=sys.stderr)
            match.annotate('EVAL', 'True')
            calls.append(match)

    calls.sort(key=lambda c: float(c.attribute('LIKESCORE')), reverse=True)
    calls = [c for c in calls if float(c.attribute('LIKESCORE')) > 0.0]
    return calls
