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

main()
{
    echo -n "Checking for wgsim: "
    which wgsim

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
}

main
