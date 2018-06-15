#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
from collections import defaultdict
import sys

import intervaltree
from intervaltree import IntervalTree
import pandas
from evalutils import IntervalForest, populate_index_from_simulation, compact
from evalutils import assess_variants_vcf, assess_variants_mvf
from evalutils import subset_variants, subset_vcf, subset_mvf
from evalutils import load_kevlar_vcf, load_triodenovo_vcf, load_gatk_mvf
import kevlar
from kevlar.vcf import VCFReader


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tolerance', type=int, metavar='T', default=10,
                        help='extend real variants by T nucleotides when '
                        'querying for overlap with variant calls; default is '
                        '10')
    parser.add_argument('--mode', choices=('Kevlar', 'GATK', 'TrioDenovo'),
                        default='Kevlar', help='Gevlar|GATK|TrioDenovo')
    parser.add_argument('--cov', default='30', help='coverage')
    parser.add_argument('--correct', help='print correct variants to file')
    parser.add_argument('--missing', help='print missing variants to file')
    parser.add_argument('--false', help='print false variants to file')
    parser.add_argument('--collisions', help='print calls that match the '
                        'same variant')
    parser.add_argument('--vartype', choices=('SNV', 'INDEL'), default=None)
    parser.add_argument('--minlength', type=int, default=None)
    parser.add_argument('--maxlength', type=int, default=None)
    parser.add_argument('--do-all', action='store_true', help='ignore all '
                        'other arguments and analyze all data')
    parser.add_argument('simvar', help='simulated variants (in custom 3-4 '
                        'column tabular format)')
    parser.add_argument('varcalls', help='VCF file of variant calls')
    return parser


def load_index(simvarfile, vartype=None, minlength=None, maxlength=None):
    with kevlar.open(simvarfile, 'r') as instream:
        if vartype:
            instream = subset_variants(
                instream, vartype, minlength=minlength, maxlength=maxlength
            )
        index = populate_index_from_simulation(instream, 'chr17')
    return index


def handle_collisions(mapping, outfile):
    numcollisions = 0
    for variant, calllist in mapping.items():
        if len(calllist) > 1:
            numcollisions += 1
    if numcollisions > 0:
        print('WARNING:', numcollisions, 'variants matched by multiple calls',
              file=sys.stderr)
        if outfile is None:
            return
        with open(outfile, 'w') as outstream:
            for variant, calllist in mapping.items():
                if len(calllist) > 1:
                    print('\n#VARIANT:', variant, file=outstream)
                    for varcall in calllist:
                        if args.mvf:
                            print('    -', varcall, file=outstream)
                        else:
                            print('    -', varcall.vcf, file=outstream)


def handle_missing(missing, outfile):
    if outfile is None:
        return
    with kevlar.open(outfile, 'w') as outstream:
        for variant in missing:
            print(variant.begin, *variant.data.split('<-'), sep='\t',
                  file=outstream)


def handle_calls(calls, outfile, mvf=False):
    if outfile is None:
        return
    with kevlar.open(outfile, 'w') as outstream:
        if mvf:
            for varcall in calls:
                print(varcall, file=outstream)
        else:
            writer = kevlar.vcf.VCFWriter(outstream)
            for varcall in calls:
                writer.write(varcall)


def evaluate(simvarfile, varcalls, mode, vartype=None, minlength=None,
             maxlength=None, tolerance=10, coverage='30', correctfile=None,
             falsefile=None, missingfile=None, collisionsfile=None):
    assert mode in ('Kevlar', 'GATK', 'TrioDenovo')
    index = load_index(simvarfile, vartype, minlength, maxlength)
    if mode == 'GATK':
        variants = load_gatk_mvf(varcalls, vartype, minlength, maxlength)
        assess_func = assess_variants_mvf
    elif mode == 'Kevlar':
        variants = load_kevlar_vcf(
            varcalls, index, delta=tolerance, vartype=vartype,
            minlength=minlength, maxlength=maxlength
        )
        assess_func = assess_variants_vcf
    elif mode == 'TrioDenovo':
        variants = load_triodenovo_vcf(
            varcalls, vartype, minlength, maxlength, coverage
        )
        assess_func = assess_variants_vcf
    correct, false, missing, mapping = assess_func(
        variants, index, delta=tolerance
    )
    
    handle_collisions(mapping, collisionsfile)
    handle_missing(missing, missingfile)
    handle_calls(correct, correctfile, mvf=(mode == 'GATK'))
    handle_calls(false, falsefile, mvf=(mode == 'GATK'))
    
    return len(mapping), len(false), len(missing)


################################################################################
def vartypestr(vartype, minlength, maxlength):
    if vartype is None:
        return 'All'
    assert vartype in ('SNV', 'INDEL')
    if vartype == 'SNV':
        return 'SNV'
    return 'INDEL {}-{}bp'.format(minlength, maxlength)


def main(args):
    correct, false, missing = evaluate(
        args.simvar, args.varcalls, args.mode, vartype=args.vartype,
        minlength=args.minlength, maxlength=args.maxlength,
        tolerance=args.tolerance, coverage=args.cov, correctfile=args.correct,
        falsefile=args.false, missingfile=args.missing,
        collisionsfile=args.collisions
    )
    vartype = vartypestr(args.vartype, args.minlength, args.maxlength)
    
    colnames = ['Caller', 'Coverage', 'VarType', 'Correct', 'False', 'Missing']
    data = [args.mode, args.cov, vartype, correct, false, missing]
    row = {c: v for c, v in zip(colnames, data)}
    table = pandas.DataFrame(columns=colnames)
    table = table.append(row, ignore_index=True)
    print(table.to_string(index=False))


def do_all():
    infiles = {
        'Kevlar': 'kevlar_calls_{cov}x.vcf.gz',
        'GATK': 'GATK_calls_{cov}x.mvf.gz',
        'TrioDenovo': 'triodenovo_calls_{cov}x.vcf.gz',
    }
    vartypes = (
        ('SNV', None, None),
        ('INDEL', 1, 10),
        ('INDEL', 11, 100),
        ('INDEL', 101, 200),
        ('INDEL', 201, 300),
        ('INDEL', 301, 400),
    )
    colnames = ['Caller', 'Coverage', 'VarType', 'Correct', 'False', 'Missing']
    table = pandas.DataFrame(columns=colnames)
    for coverage in ('10', '20', '30', '50'):
        for vartype, minlen, maxlen in vartypes:
            varstr = vartypestr(vartype, minlen, maxlen)
            for caller in ('Kevlar', 'GATK', 'TrioDenovo'):
                simvar = 'SimulatedVariants_chr17_hg38.tsv.gz'
                varcalls = infiles[caller].format(cov=coverage)
                correct, false, missing = evaluate(
                    simvar, varcalls, caller, vartype=vartype, minlength=minlen,
                    maxlength=maxlen, tolerance=args.tolerance,
                    coverage=coverage
                )
                data = [caller, coverage, varstr, correct, false, missing]
                row = {c: v for c, v in zip(colnames, data)}
                table = table.append(row, ignore_index=True)
    print(table.to_string(index=False))


if __name__ == '__main__':
    args = get_parser().parse_args()
    if args.do_all:
        do_all()
    else:
        main(args)


