#!/usr/bin/env python

from __future__ import print_function
import argparse
import re
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def printlog(debug, *args):
    if debug:
        print(*args, file=sys.stderr)


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--discard', metavar='RE', type=str, nargs='+',
                        default=None, help='regex(es) of sequence ID(s) to '
                        'discard')
    parser.add_argument('--debug', action='store_true',
                        help='print debugging output')
    parser.add_argument('--minlength', type=int, default=10000,
                        help='minimum length of sequences to keep after '
                        'removing ambiguous nucleotides; default is 10000')
    parser.add_argument('--outfile', metavar='FILE', default=sys.stdout,
                        type=argparse.FileType('w'), help='output fasta file '
                        '(prints to terminal by default)')
    parser.add_argument('infile', metavar='FILE', type=argparse.FileType('r'),
                        help='input fasta file (- for standard input)')
    return parser


def main(args):
    for record in SeqIO.parse(args.infile, 'fasta'):
        if args.discard:
            if sum([1 for rx in args.discard if re.match(rx, record.id)]) > 0:
                continue

        subseqcounter = 0
        printlog(args.debug, "DEBUG: convert to upper case", record.id)
        sequence = str(record.seq).upper()
        printlog(args.debug, "DEBUG: split seq by Ns", record.id)
        subseqs = [ss for ss in re.split('[^ACGT]+', sequence) if len(ss) > args.minlength]
        printlog(args.debug, "DEBUG: print subseqs", record.id)
        for subseq in subseqs:
            subseqcounter += 1
            subid = '{:s}_chunk_{:d}'.format(record.id, subseqcounter)
            subrecord = SeqRecord(Seq(subseq), subid, '', '')
            SeqIO.write(subrecord, args.outfile, 'fasta')


if __name__ == '__main__':
    main(get_parser().parse_args())
