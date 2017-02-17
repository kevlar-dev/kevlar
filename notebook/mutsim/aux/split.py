#!/usr/bin/env python

from __future__ import print_function
import argparse
import re
import sys

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--numfiles', type=int, metavar='N', default=1,
                    help='number of files to create')
parser.add_argument('--out', metavar='PFX', default='out',
                    help='prefix for output files (default: "out")')
parser.add_argument('infile', type=argparse.FileType('r'),
                    help='input fasta file (- for standard input)')
args = parser.parse_args()

outfiles = [open('{:s}.{:d}'.format(args.out, i + 1), 'w') for i in range(args.numfiles)]

for i, record in enumerate(SeqIO.parse(args.infile, 'fasta')):
    outfile = outfiles[i % len(outfiles)]
    SeqIO.write(record, outfile, 'fasta')
