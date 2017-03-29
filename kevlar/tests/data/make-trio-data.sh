#!/usr/bin/env bash
set -eo pipefail

simulate_reads()
{
    local refr=$1
    local sample=$2
    local seed=$3
    local trio=$4
    local erate=${5:-0.0}
    local dist=${6:-300}
    local stdev=${7:-25}
    local nreads=${8:-6000}
    local readlen=${9:-50}

    wgsim \
        -e $erate -r 0.0 -d $dist -s $stdev -N $nreads \
        -1 $readlen -2 $readlen -S $seed \
        $refr ${sample}-1.fq ${sample}-2.fq

    paste <(paste - - - - < ${sample}-1.fq) \
          <(paste - - - - < ${sample}-2.fq) \
        | tr '\t' '\n' \
        > ${trio}/${sample}.fq

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

    mkdir -p trio1/ trio2/

    # Error-free
    simulate_reads bogus-genome/refr.fa ctrl1 1234 trio1
    simulate_reads bogus-genome/refr.fa ctrl2 5678 trio1
    simulate_reads bogus-genome/seq-pool-1snp.fa case1 2468 trio1
    simulate_reads bogus-genome/seq-pool-3snps.fa case2 95616 trio1
    simulate_reads bogus-genome/seq-pool-1snp-contam.fa case3 2468 trio1
    simulate_reads bogus-genome/seq-pool-1indel.fa case4 192837 trio1
    simulate_reads chr19.order2.rand.100k.fa ctrl1 4000 trio2 0.0 570 75 11500 126
    simulate_reads chr19.order2.rand.100k.fa ctrl2 401 trio2 0.0 570 75 11500 126
    simulate_reads chr19.indel.fa case1 42 trio2 0.0 570 75 11500 126
    gzip trio2/*.fq

    # 1% error rate
    simulate_reads bogus-genome/refr.fa ctrl3 97531 trio1 0.01
    simulate_reads bogus-genome/refr.fa ctrl4 86420 trio1 0.01
    simulate_reads bogus-genome/seq-pool-1indel.fa case5 969696 trio1 0.01
    trim_reads ctrl3
    trim_reads ctrl4
    trim_reads case5

    # 2% error rate
    simulate_reads bogus-genome/refr.fa ctrl5 11111 trio1 0.02
    simulate_reads bogus-genome/refr.fa ctrl6 22222 trio1 0.02
    simulate_reads bogus-genome/seq-pool-1indel.fa case6 12121 trio1 0.02
    simulate_reads bogus-genome/seq-pool-1indel.fa case6b 98989 trio1 0.02
    trim_reads ctrl5
    trim_reads ctrl6
    trim_reads case6b

    # 3% error rate
    simulate_reads bogus-genome/refr.fa ctrl7 33333 trio1 0.03
    simulate_reads bogus-genome/refr.fa ctrl8 44444 trio1 0.03
    simulate_reads bogus-genome/seq-pool-1indel.fa case7 19191 trio1 0.03
    trim_reads ctrl7
    trim_reads ctrl8
    trim_reads case7
}

main
