#!/usr/bin/env bash
set -eo pipefail

# Check software prereqs
which wgsim
which interleave-reads.py
which trim-low-abund.py
which load-into-counting.py

# Simulate reads
wgsim -e 0.0 -r 0.0 -d 300 -s 25 -N 6000 -1 50 -2 50 -S 1234 \
    bogus-genome-refr.fa trio1-ctrl1-1.fq trio1-ctrl1-2.fq
wgsim -e 0.0 -r 0.0 -d 300 -s 25 -N 6000 -1 50 -2 50 -S 5678 \
    bogus-genome-refr.fa trio1-ctrl2-1.fq trio1-ctrl2-2.fq
wgsim -e 0.0 -r 0.0 -d 300 -s 25 -N 6000 -1 50 -2 50 -S 2468 \
    bogus-genome-seq-pool-1snp.fa trio1-case1-1.fq trio1-case1-2.fq
wgsim -e 0.0 -r 0.0 -d 300 -s 25 -N 6000 -1 50 -2 50 -S 95616 \
    bogus-genome-seq-pool-3snps.fa trio1-case2-1.fq trio1-case2-2.fq
wgsim -e 0.0 -r 0.0 -d 300 -s 25 -N 6000 -1 50 -2 50 -S 2468 \
    bogus-genome-seq-pool-1snp-contam.fa trio1-case3-1.fq trio1-case3-2.fq


# Interleave split reads, trim, and build countgraphs
for sample in "case1" "case2" "case3" "ctrl1" "ctrl2"
do
    interleave-reads.py -o - trio1-${sample}-1.fq trio1-${sample}-2.fq \
        | trim-low-abund.py -M 5e6 -k 13 -Z 10 -C 2 -V -o - - \
        | load-into-counting.py --ksize 13 -M 3e5 trio1-${sample}.counts -
done

# Interleave case reads for `kevlar find` input
for sample in "case1" "case2" "case3"
do
    interleave-reads.py -o trio1-${sample}.fq trio1-${sample}-1.fq trio1-${sample}-2.fq
done

# Cleanup
rm trio1-*-?.fq *.counts.info
