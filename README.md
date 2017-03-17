[![kevlar build status](https://img.shields.io/travis/dib-lab/kevlar.svg)](https://travis-ci.org/dib-lab/kevlar)
[![Test coverage](https://img.shields.io/codecov/c/github/dib-lab/kevlar.svg)](https://codecov.io/github/dib-lab/kevlar)
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/dib-lab/kevlar/blob/master/LICENSE)

# kevlar

---::[Reference-free variant discovery in large genomes]::---

Currently, kevlar development is focused almost entirely on finding novel variants in family-based trio and quad experimental designs.
However, the method lends itself easily to more general experimental designs, which will get more attention and support in the near future.

Although a reference genome is not required, it can be utilized to reduce data volume at an early stage in the workflow and reduce the computational demands of subsequent steps.

## Quick start

```
# Discard reads that match the reference genome perfectly --> trivially uninteresting
kevlar dump --out mother.fq mother.bam reference.fasta.gz
kevlar dump --out father.fq father.bam reference.fasta.gz
kevlar dump --out proband.fq proband.bam reference.fasta.gz

# Find "interesting" (putatively novel) k-mers in the proband
# 8G memory per sample x 3 samples = 24G
kevlar find --controls mother.fq father.fq --ctrl_max 1 --case_min 8 --ksize 27 --memory 8G --out interesting-kmers.txt proband.fq

# Filter "interesting" k-mers, group, and assemble contigs representing variants
kevlar collect --memory 1G --minabund 8 --refr reference.fasta.gz --refr-memory 4G --ksize 27 --collapse --out variants.tsv interestng-kmers.txt
```
