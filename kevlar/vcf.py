#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from datetime import date
import re
import sys
import khmer
import kevlar


class VariantAnnotationError(ValueError):
    pass


class VCFWriter(object):
    def __init__(self, source='kevlar', refr=None):
        self.samples = list()
        self._header_written = False
        self._source = source
        self._refr = refr

    def write_header(fh=sys.stdout):
        print('##fileformat=VCFv4.2', file=fh)
        print('##fileDate', date().isoformat(), sep='=', file=fh)
        print('##source', self._source, sep='=', file=fh)
        if self._refr:
            print('##reference', self._refr, sep='=', file=fh)
        print(
            '##INFO=<ID=RW,Number=1,Type=String,Description="Reference '
            'window, bounding all k-mers that contain the reference allele">',
            file=fh
        )
        print(
            '##INFO=<ID=CS,Number=1,Type=String,Description="Contig sequence '
            'used to make variant call">', file=fh
        )
        print(
            '##INFO=<ID=NC,Number=1,Type=String,Description="Explanation of '
            'why a contig alignment results in a no-call">', file=fh
        )
        print(
            '##INFO=<ID=VW,Number=1,Type=String,Description="Variant window, '
            'bounding all k-mers that contain the alternate allele">', file=fh
        )
        print(
            '##INFO=<ID=DN,Number=1,Type=Float,Description="Log likelihood '
            'that a variant is de novo">', file=fh
        )
        print(
            '##INFO=<ID=FP,Number=1,Type=Float,Description="Log likelihood '
            'that a variant is a false positive">', file=fh
        )
        print(
            '##INFO=<ID=IN,Number=1,Type=Float,Description="Log likelihood '
            'that a variant is inherited">', file=fh
        )
        print(
            '##INFO=<ID=LR,Number=1,Type=Float,Description="Likelihood ratio: '
            'DN / (FP + IN)">', file=fh
        )
        print(
            '##FORMAT=<ID=AA,Number=.,Type=Integer,Description="Abundance of '
            'all k-mers containing the alternate allele that are not present '
            'elsewhere in the genome">', file=fh
        )
        print(
            '##FORMAT=<ID=RA,Number=.,Type=Integer,Description="Abundance of '
            'reference k-mers associated with each novel alternate k-mer">',
            file=fh
        )
        print('#', end='', file=fh)
        fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        if len(self.samples) > 0:
            fields += ['FORMAT'] + self.samples
        print(*fields, sep='\t', file=fh)

    def register_sample(self, label):
        self.samples.append(label)

    def write(self, var, fh=sys.stdout):
        if not self._header_written:
            self.write_header(fh)
            self._header_written = True

        mandatory_fields = var.vcf
        format_fields = list()
        outfmt = None
        for sample in self.samples:
            fmt = list()
            values = list()
            for field in ['AA', 'RA']:
                value = var.format(sample, field)
                if value:
                    fmt.append(field)
                    values.append(value)
            fmtstr = ':'.join(fmt)
            if outfmt is None:
                outfmt = fmtstr
            else:
                if outfmt != fmtstr:
                    msg = 'samples not annotated with the same FORMAT fields'
                    msg += ' ({:s} vs {:s})'.format(outfmt, fmtstr)
                    raise VariantAnnotationError(msg)
            format_fields.append(':'.join(values))
        print(mandatory_fields, outfmt, *format_fields, sep='\t', file=fh)


class Variant(object):
    """Base class for handling variant calls and no-calls."""

    def __init__(self, seqid, pos, refr, alt, **kwargs):
        """
        Constructor method.

        The `pos` parameter expects the genomic position as a 0-based index.
        Setting the `refr` or `alt` parameters to `.` will designate this
        variant as a "no call".
        """
        self._seqid = seqid
        self._pos = pos
        self._refr = refr
        self._alt = alt
        self.info = dict()
        for key, value in kwargs.items():
            self.info[key] = value
        self.sampledata = dict()

    def format(self, sample, key, value_to_store=None):
        if value_to_store is None:
            if sample not in self.sampledata:
                return None
            if key not in self.sampledata[sample]:
                return None
            return self.sampledata[sample][key]
        else:
            self.sampledata[sample][key] = value_to_store

    def __str__(self):
        pos = self._pos + 1
        if len(self._refr) > len(self._alt):
            # Deletion
            dellength = len(self._refr) - len(self._alt)
            return '{:s}:{:d}:{:d}D'.format(self._seqid, pos, dellength)
        elif len(self._refr) < len(self._alt):
            # Insertion
            insertion = self._alt[1:]
            return '{:s}:{:d}:I->{:s}'.format(self._seqid, pos, insertion)
        else:
            # SNV
            return '{:s}:{:d}:{:s}->{:s}'.format(self._seqid, self._pos,
                                                 self._refr, self._alt)

    @property
    def seqid(self):
        return self._seqid

    @property
    def position(self):
        return self._pos

    @property
    def vcf(self):
        """Print variant to VCF."""
        attrstr = '.'
        if len(self.info) > 0:
            kvpairs = list()
            for key in sorted(self.info):
                if key != 'CS':
                    kvpairs.append(self.attribute(key, pair=True))
            queryseq = self.attribute('CS', pair=True)
            if queryseq:
                kvpairs.append(queryseq)
            attrstr = ';'.join(kvpairs)

        filterstr = 'PASS' if self._refr != '.' else '.'
        return '{:s}\t{:d}\t.\t{:s}\t{:s}\t.\t{:s}\t{:s}'.format(
            self._seqid, self._pos + 1, self._refr, self._alt, filterstr,
            attrstr
        )

    @property
    def cigar(self):
        return self.attribute('CIGAR')

    @property
    def window(self):
        """
        Getter method for the variant window.

        The "variant window" (abbreviated `VW` in VCF output) is the sequence
        interval in the proband contig that encompasses all k-mers overlapping
        the variant.

        GCCTAGTTAGCTAACGTCCCGATCACTGTGTCACTGC
                    .....A
                     ....A.
                      ...A..
                       ..A...
                        .A....
                         A.....
                         |        <-- position of variant
                    [---------]   <-- variant window, interval (inclusive)
                                      encompassing all 6-mers that overlap the
                                      variant
        """
        return self.attribute('VW')

    @property
    def refrwindow(self):
        """Similar to `window`, but encapsulating the reference allele."""
        return self.attribute('RW')

    def attribute(self, key, pair=False):
        if key not in self.info:
            return None
        value = self.info[key].replace(';', ':')
        if pair:
            keyvaluepair = '{:s}={:s}'.format(key, value)
            return keyvaluepair
        else:
            return value

    @property
    def genotypes(self):
        gt = self.attribute('GT')
        if not gt:
            return None
        return tuple(gt.split(','))
