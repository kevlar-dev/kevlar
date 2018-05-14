#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
from datetime import date
from enum import Enum
import sys
import kevlar


class VariantAnnotationError(ValueError):
    pass


class VariantFilter(Enum):
    PerfectMatch = 1
    InscrutableCigar = 2
    PassengerVariant = 3
    MateFail = 4
    PartitionScore = 5


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
        self._filters = set()
        self.info = dict()
        for key, value in kwargs.items():
            self.annotate(key, value)
        self._sample_data = defaultdict(dict)

    def __str__(self):
        if len(self._refr) == 1 and len(self._alt) == 1:
            return '{:s}:{:d}:{:s}->{:s}'.format(self._seqid, self._pos,
                                                 self._refr, self._alt)
        else:
            pos = self._pos + 1
            if len(self._refr) > len(self._alt):
                dellength = len(self._refr) - len(self._alt)
                return '{:s}:{:d}:{:d}D'.format(self._seqid, pos, dellength)
            else:
                insertion = self._alt[1:]
                return '{:s}:{:d}:I->{:s}'.format(self._seqid, pos, insertion)

    def format(self, sample, key, value_to_store=None):
        if value_to_store is None:
            if sample not in self._sample_data:
                return None
            if key not in self._sample_data[sample]:
                return None
            return self._sample_data[sample][key]
        else:
            self._sample_data[sample][key] = value_to_store

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
                if key != 'CONTIG':
                    kvpairs.append(self.attribute(key, pair=True))
            queryseq = self.attribute('CONTIG', pair=True)
            if queryseq:
                kvpairs.append(queryseq)
            attrstr = ';'.join(kvpairs)

        return '{:s}\t{:d}\t.\t{:s}\t{:s}\t.\t{:s}\t{:s}'.format(
            self._seqid, self._pos + 1, self._refr, self._alt, self.filterstr,
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
        return self.attribute('ALTWINDOW')

    @property
    def windowlength(self):
        window = self.window
        if window is None:
            return 0
        return len(window)

    @property
    def refrwindow(self):
        """Similar to `window`, but encapsulating the reference allele."""
        return self.attribute('REFRWINDOW')

    def annotate(self, key, value):
        value = str(value)
        if key in self.info:
            if isinstance(self.info[key], set):
                self.info[key].add(value)
            else:
                oldvalue = self.info[key]
                self.info[key] = set((oldvalue, value))
        else:
            self.info[key] = value

    def attribute(self, key, pair=False):
        if key not in self.info:
            return None
        value = self.info[key]
        if isinstance(value, set):
            value = ','.join(sorted(value))
        value = value.replace(';', ':')
        if pair:
            keyvaluepair = '{:s}={:s}'.format(key, value)
            return keyvaluepair
        else:
            return value

    def filter(self, filtertype):
        if not isinstance(filtertype, VariantFilter):
            return
        self._filters.add(filtertype)

    @property
    def filterstr(self):
        if len(self._filters) > 0:
            return ';'.join(sorted([vf.name for vf in self._filters]))
        elif self._refr == '.':
            return '.'
        else:
            return 'PASS'

    @property
    def genotypes(self):
        gt = self.attribute('GT')
        if not gt:
            return None
        return tuple(gt.split(','))


class VCFWriter(object):
    filter_desc = {
        VariantFilter.PerfectMatch:
            'No mismatches between contig with putatively novel content and '
            'reference target',
        VariantFilter.InscrutableCigar:
            'Alignment path/structure cannot be interpreted as a variant',
        VariantFilter.PassengerVariant:
            'A mismatch between contig and reference that is not spanned by '
            'any novel k-mers',
        VariantFilter.MateFail:
            'Aligning mate reads suggests a better location for this variant '
            'call',
        VariantFilter.PartitionScore:
            'Expectation is 1 variant call per partition, so all call(s) with '
            'suboptimal likelihood scores are filtered'
    }

    info_metadata = {
        'ALTWINDOW': (
            'String', '1', 'window containing all k-mers that span the '
            'variant alternate allele',
        ),
        'CIGAR': (
            'String', '1', 'alignment path',
        ),
        'IKMERS': (
            'Integer', '1', 'number of "interesting" (novel) k-mers spanning '
            'the variant alternate allele',
        ),
        'KSW2': (
            'Float', '1', 'alignment score',
        ),
        'MATEDIST': (
            'Float', '1', 'average distance of aligned mates of assembled '
            'novel reads',
        ),
        'REFRWINDOW': (
            'String', '1', 'window containing all k-mers that span the '
            'variant reference allele',
        ),
        'CONTIG': (
            'String', '1', 'contig assembled from reads containing novel '
            'k-mers, aligned to reference to call variants',
        ),
        'LIKESCORE': (
            'Float', '1', 'likelihood score of the variant, computed as '
            '`LLDN - max(LLIH, LLFP)`',
        ),
        'LLDN': (
            'Float', '1', 'log likelihood that the variant is a de novo '
            'variant',
        ),
        'LLIH': (
            'Float', '1', 'log likelihood that the variant is an inherited '
            'variant',
        ),
        'LLFP': (
            'Float', '1', 'log likelihood that the variant is a false call',
        ),
        'DROPPED': (
            'Integer', '1', 'number of k-mers dropped from ALTWINDOW for '
            'likelihood calculations because it is present elsewhere in the '
            'genome (not novel)',
        ),
    }

    format_metadata = {
        'ALTABUND': (
            'Integer', '.', 'abundance of alternate allele k-mers',
        ),
    }

    def __init__(self, outstream, source='kevlar', refr=None):
        self._out = outstream
        self._sample_labels = list()
        self._source = source
        self._refr = refr

    def register_sample(self, label):
        self._sample_labels.append(label)

    def describe_format(self, label, datatype, datanumber, desc):
        self.format_metadata[label] = (datatype, datanumber, desc)

    def write_header(self):
        print('##fileformat=VCFv4.2', file=self._out)
        print('##fileDate', date.today().isoformat(), sep='=', file=self._out)
        print('##source', self._source, sep='=', file=self._out)
        if self._refr:
            print('##reference', self._refr, sep='=', file=self._out)
        for filt in VariantFilter:
            filtstr = '##FILTER=<ID={label},Description="{desc}">'.format(
                label=filt.name, desc=self.filter_desc[filt]
            )
            print(filtstr, file=self._out)
        for label, (itype, inumber, idesc) in self.info_metadata.items():
            infostr = '##INFO=<ID={label},'.format(label=label)
            infostr += 'Number={num},'.format(num=inumber)
            infostr += 'Type={type},'.format(type=itype)
            infostr += 'Description="{desc}">'.format(desc=idesc)
            print(infostr, file=self._out)
        for label, (itype, inumber, idesc) in self.format_metadata.items():
            fmtstr = '##FORMAT=<ID={label},'.format(label=label)
            fmtstr += 'Number={num},'.format(num=inumber)
            fmtstr += 'Type={type},'.format(type=itype)
            fmtstr += 'Description="{desc}">'.format(desc=idesc)
            print(fmtstr, file=self._out)
        print('#', end='', file=self._out)
        fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        if len(self._sample_labels) > 0:
            fields += ['FORMAT'] + self._sample_labels
        print(*fields, sep='\t', file=self._out)

    def write(self, variant):
        fmt_fields = list()
        outfmt = None
        for sample in self._sample_labels:
            fmt = list()
            values = list()
            for field in self.format_metadata.keys():
                value = variant.format(sample, field)
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
            fmt_fields.append(':'.join(values))
        print(variant.vcf, end='', file=self._out)
        if len(fmt_fields) > 0:
            print('', outfmt, *fmt_fields, sep='\t', end='', file=self._out)
        print('\n', end='', file=self._out)


class VCFReader(object):
    def __init__(self, instream):
        self._in = instream
        self._sample_labels = list()

    def __iter__(self):
        for line in self._in:
            if not line.startswith('#CHROM\t'):
                continue
            self._save_samples(line)
            break
        for line in self._in:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            seqid = fields[0]
            pos = int(fields[1]) - 1
            refr = fields[3]
            alt = fields[4]
            variant = Variant(seqid, pos, refr, alt)
            for kvp in fields[7].split(';'):
                for key, values in kvp.split('='):
                    for value in values.split(','):
                        variant.annotate(key, value)
            if len(fields) > 9:
                fmtkeys = fields[8].split(':')
                sample_data = fields[9:]
                if len(sample_data) != len(self._sample_labels):
                    message = 'sample number mismatch: ' + line
                    raise VariantAnnotationError(message)
                for label, data in zip(self._sample_labels, sample_data):
                    fmtvalues = data.split(':')
                    if len(fmtkeys) != len(fmtvalues):
                        message = 'format data mismatch: ' + line
                        raise VariantAnnotationError(message)
                    for datakey, datavalue in zip(fmtkeys, fmtvalues):
                        variant.format(label, datakey, datavalue)
            yield variant

    def _save_samples(self, line):
        fields = line.strip().split('\t')
        assert len(fields) >= 8
        if len(fields) == 8:
            return
        self._sample_labels = fields[9:]
