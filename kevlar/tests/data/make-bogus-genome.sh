#!/usr/bin/env bash
set -eo pipefail

# Check software prereqs
which make-random-genome.py
which bwa
which samtools

mkdir -p bogus-genome

# Simulate genome
make-random-genome.py --length 15000 --seed 42 --name bogus-genome-chr1 \
    > bogus-genome/chr1.fa
make-random-genome.py --length 10000 --seed 2468 --name bogus-genome-chr2 \
    > bogus-genome/chr2.fa
make-random-genome.py --length 12000 --seed 13579 --name bogus-genome-chr3 \
    > bogus-genome/chr3.fa
cat bogus-genome/chr?.fa > bogus-genome/refr.fa

# Create genome masks
cat bogus-genome/chr2.fa bogus-genome/chr3.fa > bogus-genome/mask-chr1.fa
cat bogus-genome/chr1.fa bogus-genome/chr3.fa > bogus-genome/mask-chr2.fa

# Throw in a little contamination
make-random-genome.py --length 500 --seed 112358 --name bogus-genome-contam \
    > bogus-genome/contam1.fa
cat bogus-genome/chr?.fa bogus-genome/contam1.fa > bogus-genome/refr-contam.fa

# Generate some mutations
./mutate1.py bogus-genome/refr.fa > bogus-genome/seq-pool-1snp.fa
./mutate3.py bogus-genome/refr.fa > bogus-genome/seq-pool-3snps.fa
./mutate1.py bogus-genome/refr-contam.fa > bogus-genome/seq-pool-1snp-contam.fa
./mutate-indel.py bogus-genome/refr.fa > bogus-genome/seq-pool-1indel.fa

# Simulate and map some simple reads
./make-bogus-reads.py bogus-genome/refr-contam.fa > bogus-genome/reads-small.fq
cd bogus-genome/
bwa index refr.fa refr.fa
bwa mem -k 13 -a refr.fa reads-small.fq \
    | samtools view -b -o - \
    | samtools sort -o reads.bam -
rm refr.fa.* reads-small.fq
cd ../

# Simulate and map some reads with indels
./make-indel-reads.py bogus-genome/refr-contam.fa > bogus-genome/reads-indels.fq
cd bogus-genome/
bwa index refr.fa refr.fa
bwa mem -k 13 -a refr.fa reads-indels.fq \
    | samtools view -b -o - \
    | samtools sort -o reads-indels.bam -
rm refr.fa.* reads-indels.fq
cd ../

# Simulate an indel on random seq with same nucl composition as human chr19
./mutate-indel-chr19.py chr19.order2.rand.100k.fa > chr19.indel.fa

# Cleanup
rm -f bogus-genome/chr?.fa bogus-genome/refr-contam.fa
