*K*-mer banding
===============

If memory is a limiting factor for *k*-mer counting, kevlar supports a scatter/gather approach to based on a strategy we call "*k*-mer banding."
In brief, kevlar can achieve an N-fold reduction in memory usage in exchange for counting *k*-mers N batches.
For each batch, kevlar ignores all *k*-mers except those whose hash values fall within a specified numerical range (band), reducing the memory required to achieve accurate *k*-mer counts.

The :ref:`kevlar_count_api` and :ref:`kevlar_novel_api` commands support *k*-mer banding, and the :ref:`kevlar_unband_api` command merges novel reads from multiple batches into a single read set suitable for downstream analysis.

The example below is not intended to show the most efficient way of invoking the banding workflow, simply how it is intended to work.

.. code::

    # Count k-mers in 6 passes for a 6x reduction in memory
    for $band in {1..6}
    do
        for indiv in mother father proband
        do
            kevlar count \
                --threads 16 --memory 8G --counter-size 8 --ksize 31 \
                --mask refr-univec-mask.nodetable \
                --num-bands 6 --band $band \
                ${indiv}.kmer-counts.band${band}.counttable ${indiv}.reads.fastq.gz
        done
    done

    # Find novel k-mers in 6 passes
    for $band in {1..6}
    do
        kevlar novel \
            --case proband.reads.fastq.gz --case-counts ${indiv}.kmer-counts.band${band}.counttable \
            --control-counts mother.kmer-counts.band${band}.counttable mother.kmer-counts.band${band}.counttable \
            --case-min 5 --ctrl-max 1 --ksize 31\
            --num-bands 6 --band $band \
            --out proband-novel-${band}.augfastq.gz
    done

    # Combine the novel k-mers from all 6 passes
    kevlar unband --out proband-novel.augfastq.gz proband-novel-{1..6}.augfastq.gz

    # The rest of the kelvar simplex workflow:
    #    - kevlar filter
    #    - kevlar partition
    #    - kevlar kevlar assemble
    #    - kevlar localize
    #    - kevlar call
    #    - kevlar simlike
