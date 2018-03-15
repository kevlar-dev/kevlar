Tutorial
========

This tutorial will cover material similar to that introduced in the :doc:`kevlar quick start <quick-start>`, but with a bit more detail and commentary.


Test data
---------

The data for this tutorial comes from `a simulated 25 Mb genome <https://github.com/standage/noble>`__ and can be downloaded anonymously from the `Open Science Framework <https://osf.io/anr56/>`__.

.. highlight:: none

.. code::

    curl -L https://s3-us-west-1.amazonaws.com/noble-trios/neon-mother-reads.fq.gz -o mother.fq.gz
    curl -L https://s3-us-west-1.amazonaws.com/noble-trios/neon-father-reads.fq.gz -o father.fq.gz
    curl -L https://s3-us-west-1.amazonaws.com/noble-trios/neon-proband-reads.fq.gz -o proband.fq.gz
    curl -L https://s3-us-west-1.amazonaws.com/noble-trios/neon-refr.fa.gz -o refr.fa.gz

The reference genome must be indexed for BWA searches before proceeding.

.. code::

    bwa index refr.fa.gz


The ``kevlar simplex`` command
------------------------------

The ``kevlar simplex`` command executes the entire simplex processing workflow for discovering novel germline mutations.
The required inputs are:

- a single proband or case sample
- one or more control samples (such as parents and siblings)
- a reference genome sequence (indexed for BWA searches)
- a set of mask sequences, including the reference genome as well as any contaminants that might be present in the proband sample

This test data set doesn't include any simulated contamination, so we are simply using the reference genome sequence as the mask.

.. code::

    kevlar simplex \
         --case proband.fq.gz --control mother.fq.gz --control father.fq.gz \
         --case-min 6 --ctrl-max 1 \
         --novel-memory 500M --filter-memory 10M --filter-fpr 0.005 --mask-memory 50M  \
         --mask-files refr.fa.gz \
         --threads 4 --ksize 31 \
         --out kevlar-variant-calls.vcf \
         refr.fa.gz

Here is a discussion of various parameter choices.

- The ``--case-min`` and ``--ctrl-max`` parameters define the threshold model for selecting "interesting" (putatively novel) *k*-mers.
  Here, any *k*-mer that appears at least 6 times in the proband and at most 1 time in either parent is marked as a putative signature of a *de novo* variant.
- The ``--novel-memory`` parameter determines how much memory is allocated for computing initial *k*-mer counts in each sample.
  The accuracy of the counts depends on the number of distinct *k*-mers in the sample and the amount of memory used for counting *k*-mers.
  For human data sets sampled to about 30x coverage, allocating 10-20 Gb of memory per sample will result in high false positive rates for the initial *k*-mer counts.
  However, tests have shown that FPRs of up to 0.8 for the initial *k*-mer counts can be compensated for by subsequent filtering steps.
- The ``--filter-memory`` parameter indicates the amount of memory allocated for recomputing *k*-mer abundances of "interesting" reads.
  For human data sets sampled to about 30x coverage, 4 Gb of memory for this second filtering pass is usually more than enough for near-perfect accuracy.
  The ``--filter-fpr`` parameter indicates the level of desired accuracy (the program will terminate if the FPR is higher).


Assessing accuracy
------------------

The output of the command above is a VCF file, and in this case should contain 10 variant calls.
The ``kevlar-eval.sh`` script below uses ``bedtools intersect`` to do a quick and simple evaluation of kevlar's accuracy.

.. code:: bash

    curl -L https://s3-us-west-1.amazonaws.com/noble-trios/neon.vcf -o neon-refr.vcf
    curl -L curl -L https://raw.githubusercontent.com/standage/noble/master/kevlar-eval.sh -o kevlar-eval.sh
    bash kevlar-eval.sh neon-refr.vcf kevlar-variant-calls.vcf


The banding workflow
--------------------

If memory is a limiting factor, it's possible to reduce the most memory intensive steps of the kevlar workflow using "*k*-mer banding".
To summarize, the *k*-mer counts are computed in N independent passes over the data, each pass only requiring 1/N of the memory required for a single pass.
Banding is not yet supported in the ``kevlar simplex`` command, so the banding workflow is a bit more involved.

.. code::

    # Let's count k-mers and find novel k-mers in 6 passes for a 6x reduction in memory
    for $band in {1..6}
    do
        kevlar novel \
            --case proband.fq.gz --control mother.fq.gz --control father.fq.gz \
            --case-min 6 --ctrl-max 1 \
            --memory 500M --threads 4 --ksize 31 \
            --num-bands 6 --band $band \
            --out proband-novel-${band}.augfastq.gz
    done

    # The "kevlar filter" command will combine the results from all 6 passes and recompute k-mer counts
    kevlar filter \
        --mask refr.fa.gz --mask-memory 50M --mask-max-fpr 0.005 \
        --abund-memory 10M --abund-max-fpr 0.005 \
        --ksize 31 --out proband-filtered.augfastq.gz \
        proband-novel-{1..6}.augfastq.gz

    # The "kevlar partition" command separates the reads into sets corresponding to distinct variants.
    kevlar partition proband-partitioned.augfastq.gz proband-filtered.augfastq.gz

    # The "kevlar alac" command does assembly, alignment, and variant calling for each partition
    kevlar alac --ksize 31 --seed-size 51 --delta 50 --out proband-calls.vcf proband-partitioned.augfastq.gz refr.fa.gz


Pre-processing
--------------

Several types of pre-processing can yield large improvements in kevlar's performance.

- **Error correction**:
  The amount of memory kevlar needs for accurate variant discovery depends on the number of distinct *k*-mers in each sample, the majority of which are associated with sequencing errors.
  Using an error correction tool such as `Lighter <https://github.com/mourisl/Lighter>`__ or `BFC <https://github.com/lh3/bfc>`__ will remove many erroneous *k*-mers, reduce kevlar's memory demands without adding too much processing time to the overall workflow.
  See `this blog post <https://standage.github.io/information-content-versus-data-volume-and-k-mer-counting-accuracy.html>`__ for more details.
  However, in some cases this may lead to true variants being erroneously "corrected" and result in a loss of variant calling sensitivity.
- **Discarding reads that match the reference genome perfectly**:
  Reducing the total volume of the data set doesn't reduce the number of distinct *k*-mers to the extent that error correction does, but it can reduce the number of reads that need to be processed by up to 80%, substantially speeding up subsequent processing steps.

  .. code::

      kevlar dump --out proband-dump.fq.gz --refr refr.fa.gz proband.bam

  The time savings may not offset the time required from running ``kevlar dump`` on each sample.
  However, when testing and benchmarking kevlar on a new data set where multiple runs for parameter refinement are likely, the time savings will probably be worth the initial time investment.
- **Pre-computing count tables**:
  In the examples above, ``kevlar simplex`` and ``kevlar novel`` compute *k*-mer counts directly from Fastq data.
  However, both commands can also accept pre-computed count tables with the ``--case-counts`` and ``--control-counts`` options.
  This can be helpful when testing or benchmarking kevlar requires multiple runs over the same data.
  The ``kevlar count`` command takes Fastq files as input and saves *k*-mer count tables to disk for subsequent use.

  .. code::

      kevlar count --ksize 31 --memory 8G proband.counttable proband-paired-interleaved.fq.gz
      kevlar count --ksize 31 --memory 8G father.counttable father-r1.fq.gz father-r2.fq.gz
      kevlar count --ksize 31 --memory 8G mother.counttable mother-r1.fq.gz mother-r2.fq.gz mother-ru.fq.gz
