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


Before getting started
----------------------

For this simple example, all required input data and configuration options are provided.
However, when analyzing one's own data, it's important to consider a few points.

- **How much memory is needed to store k-mer counts for each sample?**
  Kevlar uses a probabilistic data structure to store approximate *k*-mer counts, which saves space over exact data structures but also introduces error.
  The accuracy of the approximate counts depends on two factors: the amount of memory used to store the counts, and the number of distinct *k*-mers present in the sample.
  This accuracy is summarized by a statistic called the *false positive rate (FPR)*.
  The best way to reduce Kevlar's memory requirements while keeping the FPR steady is to reduce the number of *k*-mers that need to be tracked.
  One way this is done is by creating a **mask** of uninteresting *k*-mers for kevlar to ignore: see below.
  Another way is to perform error correction on reads prior to *k*-mer counting.
  The majority of distinct *k*-mers in a sample span a sequencing error and occur only once.
  Performing error correction on the reads will eliminate many of these errors and drastically reduce the number of distinct *k*-mers present.
  Note, however, that error correction algorithms sometimes struggle to distinguish sequencing error from true variation, especially when there is low coverage.
  One can expect a small loss of sensitivity for SNV discovery if error correction is used.

  It's also important to note that different amounts of memory can be used for different samples.
  The *k*-mer counts in the case sample (proband) are recomputed and corrected throughout the workflow, and an initially high FPR can be corrected at a later stage.
  This is not true for the control samples (parents), where a high FPR will lead to a loss of sensitivity.

  We recommend an FPR ≤ 0.5 for the case sample and ≤ 0.05 for the control samples.
  For properly masked reads, this requires about **10G-20G per sample for error-corrected reads** and anywhere between **36G-72G per sample for non-error-corrected reads**.
- **Which k-mers should Kevlar mask?**
  Any *k*-mer present in the reference genome is not a signature of novel mutation, and therefore doesn't need to be tracked.
  Also, *k*-mers from any known or suspected sources of contamination (adapters, bacterial contaminants, etc.) may present false signatures of novel variation.
  Kevlar uses a **mask** to ignore such *k*-mers during counting and filtering stages of the workflow.
  It is recommended that the reference genome and any suspected sources of contamination (such as UniVec) be included in the mask.
  Consequently, if adapter sequences are in the mask, adapter trimming of reads prior to analysis is not necessary.
- **Which thresholds should I use?**
  In our experience, ``--case-min 5`` and ``--ctrl-max 1`` work well for human genomes sequenced at about 30x coverage.
  We've found that ``--ctrl-max 0`` is too strict due to sequencing errors and inflated *k*-mer counts.
  The ``--case-min`` threshold can be adjusted for higher or lower coverage, with the usual caveats: the stricter the threshold, the more confidence you have in your results but the more true variants you miss.
- **Which k-mer size should I use?**
  In our experience, *k=31* works well for both SNV and indel discovery.
  However, variants occurring in repetitive regions may be undetectable without a larger *k* size (between 45 and 55).
  Increasing the *k*-mer size will have only a moderate effect on the amount of memory required to store *k*-mer counts.
  However, note that as *k* increases each *k*-mer is more likely to span a sequencing error, and accordingly error correction will a much larger impact on performance.
- **How can I filter variant calls?**
  Kevlar implements a variety of filters to distinguish true *de novo* variants from inherited variants and false positive variant calls.
  The documentation for :ref:`kevlar_call_api` and :ref:`kevlar_simlike_api` describe these filters.
  However, it's also common to filter out variant predictions that are common variants (not *de novo*) or that occur in problematic regions such as segmental duplications or simple sequence repeats.
  The `varfilter` setting accepts a BED file that specifies genomic regions from which to filter out variant predictions.
  Make sure that the chromosome names match those in the reference genome file!


The kevlar simplex workflow: Mark I
-----------------------------------

While each step of the kevlar workflow can be invoked independently using the ``kevlar`` command, we suggest using the `Snakemake workflow <https://github.com/dib-lab/kevlar/tree/master/kevlar/workflows/mark-I>`__ provided with the kevlar source code distribution.
Instructions for configuring and running the workflow are provided `alongside the Snakefile and configuration template <https://github.com/dib-lab/kevlar/tree/master/kevlar/workflows/mark-I>`__.
The notes above should be helpful in setting parameter values and selecting appropriate input files.


Assessing accuracy
------------------

The output of the command above is a VCF file, and in this case should contain 10 variant calls.
The ``kevlar-eval.sh`` script below uses ``bedtools intersect`` to do a quick and simple evaluation of kevlar's accuracy.

.. code:: bash

    curl -L https://s3-us-west-1.amazonaws.com/noble-trios/neon.vcf -o neon-refr.vcf
    curl -L curl -L https://raw.githubusercontent.com/standage/noble/master/kevlar-eval.sh -o kevlar-eval.sh
    bash kevlar-eval.sh neon-refr.vcf kevlar-variant-calls.vcf


The kevlar simplex workflow in detail
-------------------------------------

- **[Step 0: count k-mers]** The :ref:`kevlar_count_api` command is used to count *k*-mers for each sample, as well as for the reference genome and the mask.
- **[Step 1: find interesting k-mers]** The :ref:`kevlar_novel_api` command uses pre-computed *k*-mer counts to find reads containing novel *k*-mers using the specified thresholds.
- **[Step 2: filter k-mers and reads]** The :ref:`kevlar_filter_api` command recomputes *k*-mer counts and filters out *k*-mers with insufficient abundance or *k*-mers from contaminant sources. Any reads that no longer have any interesting *k*-mers after filtering are discarded.
- **[Step 3: partition reads]** Reads spanning the same variant will typically share numerous interesting *k*-mers. The :ref:`kevlar_partition_api` command groups reads based on shared novel *k*-mers.
- **[Step 4: contig assembly]** The :ref:`kevlar_assemble_api` command assembles each partition of reads into contigs for variant annotation.
- **[Step 5: localize reference targets]** The :ref:`kevlar_localize_api` command identifies the appropriate target (or set of targets) in the reference genome for aligning each variant-spanning contig for variant annotation.
- **[Step 6: call variants]** The :ref:`kevlar_call_api` command computes a full dynamic programming alignment of each reference-spanning contig to its corresponding reference target(s) and calls variants based on the alignment path.
- **[Step 7: score and rank variant calls]** The :ref:`kevlar_simlike_api` command computes a likelihood score for each variant prediction and ranks variant calls based on this score.
