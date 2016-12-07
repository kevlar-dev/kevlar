#!/usr/bin/env bash
set -eo pipefail

# Check software prereqs
which make-random-genome.py

# Simulate genome
make-random-genome.py --length 15000 --seed 42 --name bogus-genome-chr1 \
    > bogus-genome-chr1.fa
make-random-genome.py --length 10000 --seed 2468 --name bogus-genome-chr2 \
    > bogus-genome-chr2.fa
make-random-genome.py --length 12000 --seed 13579 --name bogus-genome-chr3 \
    > bogus-genome-chr3.fa
make-random-genome.py --length 5000 --seed 2048 --name bogus-genome-contam \
    > bogus-genome-contam.fa

cat bogus-genome-chr?.fa > bogus-genome-refr.fa
cat bogus-genome-chr?.fa bogus-genome-contam.fa > bogus-genome-seq-pool.fa
./mutate1.py bogus-genome-seq-pool.fa > bogus-genome-seq-pool-1snp.fa
./mutate4.py bogus-genome-seq-pool.fa > bogus-genome-seq-pool-4snps.fa
cat bogus-genome-chr2.fa bogus-genome-chr3.fa > bogus-genome-mask-chr1.fa
cat bogus-genome-chr1.fa bogus-genome-chr3.fa > bogus-genome-mask-chr2.fa
