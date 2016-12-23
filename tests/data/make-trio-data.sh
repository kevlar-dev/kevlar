#!/usr/bin/env bash
set -eo pipefail

simulate_reads()
{
    local refr=$1
    local sample=$2
    local seed=$3
    local erate=${4:-0.0}
    local dist=${5:-300}
    local stdev=${6:-25}
    local nreads=${7:-6000}
    local readlen=${8:-50}

    wgsim \
        -e $erate -r 0.0 -d $dist -s $stdev -N $nreads \
        -1 $readlen -2 $readlen -S $seed \
        $refr ${sample}-1.fq ${sample}-2.fq
    
    paste <(paste - - - - < ${sample}-1.fq) \
          <(paste - - - - < ${sample}-2.fq) \
        | tr '\t' '\n' \
        > trio1/${sample}.fq

    rm ${sample}-[1,2].fq
}

trim_reads()
{
    local sample=$1
    trim-low-abund.py -M 5e6 -k 13 -Z 10 -C 2 -V -o ${sample}.tmp trio1/${sample}.fq 
    mv ${sample}.tmp trio1/${sample}.fq
}

main()
{
    echo -n "Checking for wgsim: "
    which wgsim
    echo -n "Checking for trim-low-abund.py: "
    which trim-low-abund.py

    mkdir -p trio1/

    # Error-free
    simulate_reads bogus-genome/refr.fa ctrl1 1234
    simulate_reads bogus-genome/refr.fa ctrl2 5678
    simulate_reads bogus-genome/seq-pool-1snp.fa case1 2468
    simulate_reads bogus-genome/seq-pool-3snps.fa case2 95616
    simulate_reads bogus-genome/seq-pool-1snp-contam.fa case3 2468
    simulate_reads bogus-genome/seq-pool-1indel.fa case4 192837

    # 1% error rate
    simulate_reads bogus-genome/refr.fa ctrl3 97531 0.01
    simulate_reads bogus-genome/refr.fa ctrl4 86420 0.01
    simulate_reads bogus-genome/seq-pool-1indel.fa case5 969696 0.01
    trim_reads ctrl3
    trim_reads ctrl4
    trim_reads case5

    # 3% error rate
    simulate_reads bogus-genome/refr.fa ctrl5 11111 0.03
    simulate_reads bogus-genome/refr.fa ctrl6 22222 0.03
    simulate_reads bogus-genome/seq-pool-1indel.fa case6 33333 0.03
    trim_reads ctrl5
    trim_reads ctrl6
    trim_reads case6
}

main
