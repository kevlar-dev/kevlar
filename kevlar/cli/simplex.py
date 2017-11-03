#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from khmer import khmer_args


def subparser(subparsers):
    desc = 'Call novel germline mutations in a simplex pedigree'
    subparser = subparsers.add_parser(
        'simplex', add_help=False, description=desc
    )
    subparser._positionals.title = 'Required inputs'
    subparser.add_argument(
        'refr', help='reference genome in Fasta format (indexed for bwa '
        'search)'
    )

    novel_args = subparser.add_argument_group(
        'Find reads containing putatively novel k-mers',
        'The "kevlar novel" procedure identifies "interesting" k-mers, or '
        'putatively novel k-mers that are high abundance in the case and '
        'absent in the controls.'
    )
    novel_args.add_argument(
        '--case', metavar='F', nargs='+', required=True, help='one or more '
        'FASTA/FASTQ files containing reads from the case/proband sample')
    novel_args.add_argument(
        '--case-counts', metavar='F', help='a counttable file containing pre-'
        'computed k-mer abundances for the case sample; if not provided, '
        'k-mer abundances will be computed from the FASTA/FASTQ input')
    novel_args.add_argument(
        '--control', metavar='F', nargs='+', action='append',
        help='one or more FASTA/FASTQ files containing reads from a control '
        'sample; must be declared once for each control sample'
    )
    novel_args.add_argument(
        '--control-counts', metavar='F', nargs='+', help='list of counttable '
        'file(s) corresponding to each control sample; if not provided, k-mer '
        'abundances will be computed from FASTA/FASTQ input'
    )
    novel_args.add_argument(
        '-x', '--ctrl-max', metavar='X', type=int, default=1,
        help='k-mers with abund > X in any control sample are uninteresting; '
        'default is X=1'
    )
    novel_args.add_argument(
        '-y', '--case-min', metavar='Y', type=int, default=5, help='k-mers '
        'with abund < Y in the case sample are uninteresting or unreliable; '
        'default is Y=5'
    )
    novel_args.add_argument(
        '--novel-memory', default='1e6', type=khmer_args.memory_setting,
        metavar='MEM', help='total memory allocated to k-mer abundance for '
        'each sample; default is 1M; ignored when pre-computed k-mer '
        'abundances are supplied via counttable'
    )
    novel_args.add_argument(
        '--novel-fpr', type=float, default=0.2, metavar='FPR',
        help='terminate if the expected false positive rate for any sample is '
        'higher than the specified FPR; default is 0.2'
    )
    novel_args.add_argument(
        '-t', '--threads', type=int, default=1, metavar='T', help='number of '
        'threads to use for processing FASTA/FASTQ files; default is 1'
    )

    filter_args = subparser.add_argument_group(
        'Filter out uninteresting k-mers or reads',
        'Some reads initially labeled as "interesting" are not in fact novel. '
        'Some contain k-mers with inflated abundance that, when corrected, do '
        'not satisfy the --case-min threshold. Other k-mers are present in '
        'the reference genome or a contaminant database and should be '
        'ignored. After the initial pass over all of the data, the abundances '
        'of the novel k-mers are recomputed with a nearly perfect accuracy '
        'with much less memory.'
    )
    filter_args.add_argument(
        '--filter-memory', default='1e6', type=khmer_args.memory_setting,
        metavar='MEM', help='memory to allocate for recomputing abundances of '
        'novel k-mers; default is 1M'
    )
    filter_args.add_argument(
        '--mask-files', metavar='MSK', type=str, default=None, nargs='+',
        help='sequences to mask (reference genomes, contaminants); can '
        'provide as one or more Fasta/Fastq files or a single pre-computed '
        'nodetable file'
    )
    filter_args.add_argument(
        '--mask-memory', default='1e6', type=khmer_args.memory_setting,
        metavar='MEM', help='memory allocated to storing k-mer abundances in '
        'the mask; default is 1M; ignored when pre-computed k-mer abundances '
        'are supplied via counttable'
    )
    filter_args.add_argument(
        '--filter-fpr', type=float, default=0.001, metavar='FPR',
        help='terminate if the expected false positive rate for the '
        'recomputed k-mer counts or the mask is higher than the specified '
        'FPR; default is 0.001'
    )

    partition_args = subparser.add_argument_group(
        'Partitioning',
        'Reads containing novel k-mers are partitioned into sets '
        'corresponding to distinct variants.'
    )
    partition_args.add_argument(
        '--part-min-abund', type=int, default=2, metavar='X', help='ignore '
        'k-mers with abundance lower than X; default is 2'
    )
    partition_args.add_argument(
        '--part-max-abund', type=int, default=200, metavar='Y', help='ignore '
        'k-mers with abundance higher than Y; default is 200'
    )
    partition_args.add_argument(
        '--no-dedup', dest='dedup', action='store_false', default=True,
        help='skip steps to remove PCR duplicates'
    )

    alac_args = subparser.add_argument_group(
        'Assembling, aligning, and calling each variant',
        'For each partition, the reads are assembled into a contig, the '
        'contig is aligned to a cutout of the reference genome, and the '
        'variant is called from the alignment.'
    )
    alac_args.add_argument(
        '-d', '--delta', type=int, default=50, metavar='D', help='extend the '
        'genomic cutout by D bp; default is 50'
    )
    alac_args.add_argument(
        '-A', '--match', type=int, default=1, metavar='A', help='alignment '
        'match score; default is 1'
    )
    alac_args.add_argument(
        '-B', '--mismatch', type=int, default=2, metavar='B', help='alignment '
        'mismatch penalty; default is 2'
    )
    alac_args.add_argument(
        '-O', '--open', type=int, default=5, metavar='O', help='alignment gap '
        'open penalty; default is 5'
    )
    alac_args.add_argument(
        '-E', '--extend', type=int, default=0, metavar='E', help='alignment '
        'gap extension penalty; default is 0'
    )

    misc_args = subparser.add_argument_group('Miscellaneous settings')
    misc_args.add_argument(
        '-h', '--help', action='help', help='show this help message and exit'
    )
    misc_args.add_argument(
        '-o', '--out', metavar='FILE', help='output file; default is terminal '
        '(stdout)'
    )
    misc_args.add_argument(
        '-k', '--ksize', type=int, default=31, metavar='K', help='k-mer size; '
        'default is 31'
    )
